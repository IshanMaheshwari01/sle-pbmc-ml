#!/usr/bin/env python
"""Phase 1, step 1: quality-control check on the cell-level subsample.

The authors already filtered cells and removed doublets, so this step checks
their work rather than starting from scratch. Thresholds adapt to each
processing batch (median +/- N MADs) instead of one fixed cut-off, because
batches were sequenced to different depths.

Doublet detection is deliberately NOT re-run here: tools such as Scrublet
estimate doublets from the full set of cells in a 10x lane, and this object
holds a 300-cells-per-sample subsample, which would distort that estimate.

Outputs
    data/processed/<subsample>_qc.h5ad        cells passing QC, raw counts
    results/reports/03_qc.md                  thresholds and what was removed
    results/figures/03_qc_by_batch.png

Usage:  python scripts/03_qc.py [--config config/config.yaml] [--dev]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.phase1 import gene_symbols, load_phase1, mad_outlier, phase1_paths, write_report  # noqa: E402
from sle_pbmc.plotting import INK_2, SLOT, apply_style, plt  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--phase1-config", default=None)
    ap.add_argument("--dev", action="store_true", help="use the tiny dev subsample")
    args = ap.parse_args()

    cfg, p1 = load_config(args.config), load_phase1(args.phase1_config)
    paths = phase1_paths(cfg, dev=args.dev)
    cols = cfg["columns"]
    batch = p1["batch_key"]
    q = p1["qc"]
    if not paths["input"].exists():
        sys.exit(f"{paths['input']} not found. Run `make subsample` (or `make dev`) first.")

    adata = ad.read_h5ad(paths["input"])
    sym = gene_symbols(adata.var, cols["gene_symbol"])
    adata.var["mt"] = sym.str.upper().str.startswith("MT-")
    adata.var["ribo"] = sym.str.upper().str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True)

    obs = adata.obs
    flags = pd.DataFrame(index=obs.index)
    flags["low_genes"] = obs["n_genes_by_counts"] < q["min_genes"]
    flags["counts_outlier"] = False
    flags["genes_outlier"] = False
    flags["mt_high"] = obs["pct_counts_mt"] > q["max_pct_mt"]
    for b, idx in obs.groupby(batch, observed=True).groups.items():
        o = obs.loc[idx]
        flags.loc[idx, "counts_outlier"] = mad_outlier(o["log1p_total_counts"], q["mad_log_counts"]).values
        flags.loc[idx, "genes_outlier"] = mad_outlier(o["log1p_n_genes_by_counts"], q["mad_log_genes"]).values
        flags.loc[idx, "mt_high"] |= mad_outlier(o["pct_counts_mt"], q["mad_pct_mt"], side="upper").values
    flags["any"] = flags.any(axis=1)
    adata.obs["qc_pass"] = ~flags["any"].values

    # ---- report -----------------------------------------------------------
    lines = ["# Phase 1 QC check\n"]
    lines.append(f"- Input: `{paths['input'].name}`, {adata.n_obs:,} cells")
    lines.append(
        f"- Rules (per batch `{batch}`): log total counts or log genes > {q['mad_log_counts']}/"
        f"{q['mad_log_genes']} MADs from the batch median; % mitochondrial > {q['mad_pct_mt']} MADs "
        f"above the median or > {q['max_pct_mt']}%; fewer than {q['min_genes']} genes."
    )
    n_mt = int(adata.var["mt"].sum())
    lines.append(f"- Mitochondrial genes found: {n_mt}" + ("" if n_mt else " (none: % mito check is inactive)"))
    removed = int(flags["any"].sum())
    lines.append(f"- **Flagged: {removed:,} cells ({removed / adata.n_obs:.2%})**; kept {adata.n_obs - removed:,}\n")

    by_batch = flags.groupby(obs[batch].astype(str).values).mean().mul(100).round(2)
    by_batch.insert(0, "cells", obs.groupby(batch, observed=True).size().values)
    lines.append("## % of cells flagged, by batch and rule\n")
    lines.append(by_batch.rename_axis("batch").to_markdown())

    med = obs.groupby(batch, observed=True)[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]].median()
    lines.append("\n## Median QC metrics by batch\n")
    lines.append(med.round(1).rename_axis("batch").to_markdown())

    by_dis = flags["any"].groupby(obs[cols["disease"]].astype(str).values).mean().mul(100).round(2)
    lines.append("\n## % flagged by disease\n")
    lines.append(by_dis.rename("% flagged").to_frame().to_markdown())
    lines.append(
        "\nIf one disease group loses many more cells than the other, QC itself would "
        "introduce a difference between groups; similar rates mean it does not."
    )

    per_sample = flags["any"].groupby(obs["sample_id"].astype(str).values).mean().mul(100)
    lines.append("\n## Samples with the highest % flagged\n")
    lines.append(per_sample.sort_values(ascending=False).head(5).round(1).rename("% flagged").to_frame().to_markdown())

    # ---- figure ----------------------------------------------------------------
    apply_style()
    metrics = [
        ("total_counts", "UMI counts per cell", True),
        ("n_genes_by_counts", "Genes detected per cell", True),
        ("pct_counts_mt", "% mitochondrial counts", False),
    ]
    batches = sorted(obs[batch].astype(str).unique())
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
    for ax, (col, label, logy) in zip(axes, metrics):
        data = [obs.loc[obs[batch].astype(str) == b, col].to_numpy() for b in batches]
        bp = ax.boxplot(
            data, widths=0.5, showfliers=False, patch_artist=True,
            medianprops=dict(color="#0b0b0b", linewidth=1.2),
            boxprops=dict(facecolor=SLOT[0], edgecolor=SLOT[0], alpha=0.85),
            whiskerprops=dict(color=INK_2), capprops=dict(color=INK_2),
        )
        ax.set_xticks(range(1, len(batches) + 1), batches)
        ax.set_xlabel(f"Processing batch ({batch})")
        ax.set_title(label)
        if logy:
            from matplotlib.ticker import FuncFormatter, LogLocator

            ax.set_yscale("log")
            ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    fig.suptitle("Cell quality by processing batch (boxes: middle 50% of cells)", fontsize=10, x=0.01, ha="left", y=1.04)
    paths["figures"].mkdir(parents=True, exist_ok=True)
    fig_path = paths["figures"] / "03_qc_by_batch.png"
    fig.savefig(fig_path)
    plt.close(fig)
    lines.append(f"\nFigure: `{fig_path.name}`")

    write_report(paths["reports"] / "03_qc.md", lines)

    if removed == adata.n_obs:
        sys.exit("QC flagged every cell; the thresholds in config/phase1.yaml do not fit this data.")
    if removed / adata.n_obs > 0.2:
        print(f"[warn] QC flagged {removed / adata.n_obs:.0%} of cells. The authors already filtered "
              "this data, so that is unexpected; check 03_qc.md before going on.")
    out = adata[adata.obs["qc_pass"].values].copy()
    out.write_h5ad(paths["qc"], compression="gzip")
    print(f"Wrote {paths['qc']}: {out.n_obs:,} cells")


if __name__ == "__main__":
    main()
