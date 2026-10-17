#!/usr/bin/env python
"""One streaming pass over the 12 GB file that produces two objects.

1. A cell-level subsample (raw counts): every sample is kept, large samples
   are thinned to `max_cells_per_sample`. Used for QC plots, integration
   (scVI/Harmony) and UMAPs, which do not need all 1.26M cells.

2. Pseudobulk counts summed per sample x cell type over ALL cells. Used for
   differential expression and the donor-level classifier, where the donor is
   the unit of replication and throwing cells away would only cost precision.

Peak memory is roughly one block of rows plus the subsample, so the default
config fits a 16 GB laptop.

Usage:
    python scripts/02_subsample.py            # full subsample (config: subsample)
    python scripts/02_subsample.py --dev      # tiny subsample for code iteration
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.h5ad_stream import iter_row_blocks, locate_counts, read_obs, read_var  # noqa: E402
from sle_pbmc.sampling import (  # noqa: E402
    PseudobulkAccumulator,
    make_sample_ids,
    select_cells,
    summarise_groups,
)


def clean_obs(obs: pd.DataFrame) -> pd.DataFrame:
    """Drop unused categories so the written h5ad does not carry empty levels."""
    obs = obs.copy()
    for col in obs.columns:
        if isinstance(obs[col].dtype, pd.CategoricalDtype):
            obs[col] = obs[col].cat.remove_unused_categories()
    return obs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--dev", action="store_true", help="tiny subsample for fast iteration")
    ap.add_argument("--skip-pseudobulk", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    cols = cfg["columns"]
    path = Path(cfg["data"]["raw_h5ad"])
    if not path.exists():
        sys.exit(f"{path} not found. Run `make download` first.")
    sub_cfg = cfg["dev"] if args.dev else cfg["subsample"]
    seed = cfg["subsample"]["seed"]
    block_size = cfg["subsample"]["block_size"]

    t0 = time.time()
    print("Reading cell metadata ...")
    obs = read_obs(path)
    loc = locate_counts(path)
    var = read_var(path, loc)
    n_cells, n_genes = loc.shape
    if len(obs) != n_cells:
        sys.exit(f"obs has {len(obs)} rows but the count matrix has {n_cells}.")

    sample_ids = make_sample_ids(obs, cfg["sample_key"])
    keep = select_cells(sample_ids, sub_cfg["max_cells_per_sample"], seed)
    keep_mask = np.zeros(n_cells, dtype=bool)
    keep_mask[keep] = True
    print(
        f"{n_cells:,} cells in {sample_ids.nunique()} samples -> keeping {keep.size:,} "
        f"(<= {sub_cfg['max_cells_per_sample']} per sample)"
    )

    acc = None
    if not args.skip_pseudobulk:
        # A pseudobulk profile is one sample x cell type, optionally split further
        # (split_by). Splitting by processing batch keeps a sample's technical
        # replicates from different batches as separate profiles, so batch stays
        # a clean covariate and those replicates can be used to estimate it.
        split_by = list(cfg["pseudobulk"].get("split_by") or [])
        unit_ids = make_sample_ids(obs, list(cfg["sample_key"]) + split_by)
        group_ids = unit_ids + "||" + obs[cols["cell_type"]].astype(str).values
        codes, uniques = pd.factorize(group_ids, sort=True)
        acc = PseudobulkAccumulator(codes, len(uniques), n_genes)

    kept_blocks = []
    n_blocks = -(-n_cells // block_size)
    for start, block in tqdm(
        iter_row_blocks(path, loc, block_size), total=n_blocks, desc="Streaming", unit="block"
    ):
        if start == 0 and block.nnz and not np.allclose(block.data, np.rint(block.data)):
            print(
                f"[warn] {loc.matrix_path} holds non-integer values; pseudobulk will round "
                "them, but DESeq2/edgeR need true raw counts. Check `make inspect`."
            )
        local = keep_mask[start : start + block.shape[0]]
        if local.any():
            kept_blocks.append(block[np.flatnonzero(local)].astype(np.float32))
        if acc is not None:
            acc.add_block(start, block)

    # ---- cell-level subsample ------------------------------------------
    X = sp.vstack(kept_blocks, format="csr")
    sub_obs = clean_obs(obs.iloc[keep])
    sub_obs["sample_id"] = pd.Categorical(sample_ids.iloc[keep].values)
    sub = ad.AnnData(X=X, obs=sub_obs, var=var.copy())
    sub.layers["counts"] = sub.X.copy()
    sub.uns["provenance"] = {
        "source_file": path.name,
        "source_cells": int(n_cells),
        "max_cells_per_sample": int(sub_cfg["max_cells_per_sample"]),
        "seed": int(seed),
        "counts_from": loc.matrix_path,
    }
    out = Path(sub_cfg["out_h5ad"])
    out.parent.mkdir(parents=True, exist_ok=True)
    sub.write_h5ad(out, compression="gzip")
    print(f"Wrote {out}: {sub.n_obs:,} cells x {sub.n_vars:,} genes")
    del X, kept_blocks, sub  # free RAM before building the pseudobulk object

    # ---- pseudobulk ------------------------------------------------------
    if acc is not None:
        sums = acc.result()
        pb_obs = summarise_groups(obs, group_ids, sample_ids, cols)
        pb_obs = pb_obs.loc[uniques]  # same order as the accumulator rows
        unit_of_group = pd.Series(unit_ids.values, index=group_ids.values)
        unit_of_group = unit_of_group[~unit_of_group.index.duplicated()]
        pb_obs.insert(1, "profile_id", unit_of_group.loc[pb_obs.index].values)
        pb_obs["low_cell_count"] = pb_obs["n_cells"] < cfg["pseudobulk"]["min_cells"]
        if sums.nnz and sums.data.max() >= np.iinfo(np.int32).max:
            sys.exit("A pseudobulk sum exceeds int32; store as int64 instead.")
        pb = ad.AnnData(
            X=sums.astype(np.int32),
            obs=pb_obs,
            var=var.copy(),
        )
        pb.obs.index = pb.obs.index.astype(str)
        pb.uns["provenance"] = {
            "source_file": path.name,
            "cells_aggregated": int(n_cells),
            "group_by": list(cfg["sample_key"]) + split_by + [cols["cell_type"]],
        }
        pb_out = Path(cfg["pseudobulk"]["out_h5ad"])
        pb_out.parent.mkdir(parents=True, exist_ok=True)
        pb.write_h5ad(pb_out, compression="gzip")
        print(
            f"Wrote {pb_out}: {pb.n_obs:,} sample x cell-type profiles "
            f"({int(pb_obs['low_cell_count'].sum())} built from < {cfg['pseudobulk']['min_cells']} cells)"
        )

    print(f"Done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
