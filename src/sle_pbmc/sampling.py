"""Sample definitions, donor-preserving subsampling and pseudobulk aggregation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def make_sample_ids(obs: pd.DataFrame, sample_key: list[str]) -> pd.Series:
    """One ID per biological sample (donor x time point by default)."""
    missing = [c for c in sample_key if c not in obs.columns]
    if missing:
        raise KeyError(f"sample_key columns missing from obs: {missing}")
    parts = [obs[c].astype(str) for c in sample_key]
    sid = parts[0]
    for p in parts[1:]:
        sid = sid + "|" + p
    return sid.rename("sample_id")


def select_cells(
    sample_ids: pd.Series, max_per_sample: int, seed: int
) -> np.ndarray:
    """Return sorted row indices: every sample kept, each capped at max_per_sample.

    Thinning is done within each sample, so no donor disappears and the
    disease/healthy balance at the donor level is unchanged. Cell-type
    proportions inside a sample are preserved in expectation (simple random
    sampling within the sample).
    """
    rng = np.random.default_rng(seed)
    codes, _ = pd.factorize(sample_ids, sort=True)
    keep = []
    for code in np.unique(codes):
        rows = np.flatnonzero(codes == code)
        if rows.size > max_per_sample:
            rows = rng.choice(rows, size=max_per_sample, replace=False)
        keep.append(rows)
    return np.sort(np.concatenate(keep))


class PseudobulkAccumulator:
    """Sum raw counts per group while streaming row blocks.

    Summed over hundreds of cells, a pseudobulk profile is close to dense, so
    the running total is a dense int64 array (n_groups x n_genes): about
    0.9 GB for ~3,900 sample x cell-type groups and 30k genes, and it never
    grows. Each block is collapsed to per-group sums (G @ X_block, at most as
    many non-zeros as the block itself) and scattered into that array.
    Integer arithmetic keeps the sums exact.
    """

    def __init__(self, group_codes: np.ndarray, n_groups: int, n_genes: int):
        self.group_codes = np.asarray(group_codes)
        self.n_groups = n_groups
        self.sums = np.zeros((n_groups, n_genes), dtype=np.int64)

    def add_block(self, start: int, block: sp.csr_matrix) -> None:
        codes = self.group_codes[start : start + block.shape[0]]
        valid = codes >= 0
        if not valid.any():
            return
        G = sp.csr_matrix(
            (np.ones(int(valid.sum()), dtype=np.int64), (codes[valid], np.flatnonzero(valid))),
            shape=(self.n_groups, block.shape[0]),
        )
        Xi = block.copy()
        Xi.data = np.rint(Xi.data)
        part = (G @ Xi.astype(np.int64)).tocoo()  # one entry per (group, gene)
        self.sums[part.row, part.col] += part.data

    def result(self) -> sp.csr_matrix:
        return sp.csr_matrix(self.sums)


def summarise_groups(
    obs: pd.DataFrame,
    group_ids: pd.Series,
    sample_ids: pd.Series,
    cols: dict[str, str],
) -> pd.DataFrame:
    """Per-group metadata for the pseudobulk object.

    For each sample-level column, a group gets the value if it is constant
    across its cells, or "mixed" otherwise (e.g. a sample whose technical
    replicates were run in two processing batches). Mixed values are reported
    by the audit, never silently picked.
    """
    df = obs.copy()
    df["_group"] = group_ids.values
    df["sample_id"] = sample_ids.values
    g = df.groupby("_group", observed=True, sort=True)

    out = pd.DataFrame(index=g.size().index)
    out.index.name = "group_id"
    out["n_cells"] = g.size()
    out["sample_id"] = g["sample_id"].first()

    def constant_or_mixed(col: str) -> pd.Series:
        nun = g[col].nunique(dropna=False)
        first = g[col].first().astype(str)
        return first.where(nun <= 1, "mixed")

    for key in ("donor", "disease", "disease_state", "batch", "ancestry", "sex", "cell_type"):
        c = cols.get(key)
        if c and c in df.columns:
            out[key] = constant_or_mixed(c)
    if cols.get("library") in df.columns:
        out["n_libraries"] = g[cols["library"]].nunique()
    if cols.get("batch") in df.columns:
        out["n_batches"] = g[cols["batch"]].nunique()
    return out
