#!/usr/bin/env python
"""Phase 1, step 3: check the authors' cell-type labels with CellTypist.

CellTypist is a logistic-regression classifier trained on a large immune
atlas. Its `Immune_All_Low` model predicts ~100 fine-grained immune types, so
besides checking the 11 author labels it also splits them further (for example
helper T cells into naive/central-memory, effector-memory and regulatory).

Each CellTypist label is mapped to a broad lineage so the two annotations can
be compared: "agreement" is the share of cells whose CellTypist lineage matches
the author label. Majority voting smooths predictions over the Leiden clusters
from the scVI graph.

Outputs
    results/reports/05_annotation.md
    results/reports/05_celltypist_predictions.csv
    results/figures/05_celltypist_vs_author.png

Usage:  python scripts/05_annotate.py [--config ...] [--dev] [--model-path my_model.pkl]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.phase1 import gene_symbols, load_phase1, phase1_paths, write_report  # noqa: E402
from sle_pbmc.plotting import INK, SURFACE, apply_style, plt  # noqa: E402

# Author labels (CELLxGENE ontology names) -> broad lineage
AUTHOR_LINEAGE = {
    "CD4-positive, alpha-beta T cell": "CD4 T",
    "CD8-positive, alpha-beta T cell": "CD8 T",
    "B cell": "B",
    "natural killer cell": "NK",
    "classical monocyte": "classical monocyte",
    "non-classical monocyte": "non-classical monocyte",
    "plasmacytoid dendritic cell": "pDC",
    "conventional dendritic cell": "cDC",
    "plasmablast": "plasmablast",
    "progenitor cell": "progenitor",
    "lymphocyte": "lymphocyte (unspecified)",
}

# blue sequential ramp (light -> dark), from the validated palette
BLUE_RAMP = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def celltypist_lineage(label: str) -> str:
    """Map a CellTypist Immune_All label to a broad lineage by keyword."""
    s = label.lower()
    if "non-classical monocyte" in s:
        return "non-classical monocyte"
    if "classical monocyte" in s:
        return "classical monocyte"
    if "monocyte" in s:
        return "monocyte (other)"
    if "helper" in s or "regulatory t" in s or "treg" in s or "follicular" in s:
        return "CD4 T"
    if "cytotoxic" in s or "mait" in s or "gamma-delta" in s or "cd8" in s:
        return "CD8 T"
    if "nk" in s or "ilc" in s:
        return "NK"
    if "plasma" in s:
        return "plasmablast"
    if "b cell" in s or "b cells" in s:
        return "B"
    if "pdc" in s:
        return "pDC"
    if "dc" in s or "dendritic" in s:
        return "cDC"
    if "hsc" in s or "mpp" in s or "progenitor" in s:
        return "progenitor"
    if "megakaryocyte" in s or "platelet" in s:
        return "platelet"
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--phase1-config", default=None)
    ap.add_argument("--dev", action="store_true")
    ap.add_argument("--model-path", default=None, help="use a local CellTypist .pkl instead of downloading")
    args = ap.parse_args()

    import celltypist
    from celltypist import models

    cfg, p1 = load_config(args.config), load_phase1(args.phase1_config)
    paths = phase1_paths(cfg, dev=args.dev)
    cols = cfg["columns"]
    if not paths["integrated"].exists():
        sys.exit(f"{paths['integrated']} not found. Run `make integrate` first.")

    adata = ad.read_h5ad(paths["integrated"])
    # CellTypist wants log1p(CP10K) expression with gene symbols as names; X already holds that.
    q = ad.AnnData(X=adata.X, obs=adata.obs[["leiden_scvi"]].copy())
    q.var_names = gene_symbols(adata.var, cols["gene_symbol"]).astype(str)
    q.var_names_make_unique()

    if args.model_path:
        model = models.Model.load(model=args.model_path)
        model_name = Path(args.model_path).name
    else:
        model_name = p1["celltypist"]["model"]
        models.download_models(model=[model_name])
        model = models.Model.load(model=model_name)

    # Predict in chunks: CellTypist densifies its input, and 80k+ cells at once
    # would need several GB. Majority voting is then done here: each Leiden
    # cluster takes its most common per-cell prediction (what CellTypist's own
    # majority_voting does with a given over-clustering).
    chunk = int(p1["celltypist"].get("chunk_size", 20000))
    labels, confs = [], []
    for start in range(0, q.n_obs, chunk):
        part = q[start : start + chunk].copy()
        pr = celltypist.annotate(part, model=model, majority_voting=False)
        r = pr.to_adata(insert_conf=True).obs
        labels.append(r["predicted_labels"].astype(str).to_numpy())
        confs.append(r["conf_score"].to_numpy() if "conf_score" in r else np.full(part.n_obs, np.nan))
        print(f"CellTypist: {min(start + chunk, q.n_obs):,}/{q.n_obs:,} cells")
    per_cell = pd.Series(np.concatenate(labels), index=q.obs_names)
    n_overclust = None
    if p1["celltypist"]["majority_voting"]:
        # Vote over a deliberately fine clustering, as CellTypist does: with a
        # coarse one, rare populations (plasmablasts, progenitors) are swallowed
        # by a large neighbouring cluster and inherit its label.
        import scanpy as sc

        res = p1["celltypist"].get("overcluster_resolution", "auto")
        if res == "auto":  # CellTypist's own rule of thumb
            n = adata.n_obs
            res = 5 if n < 5000 else 10 if n < 20000 else 15 if n < 40000 else 20 if n < 100000 else 25
        sc.tl.leiden(adata, neighbors_key="scVI", resolution=float(res), key_added="overclust",
                     flavor="igraph", n_iterations=2, directed=False, random_state=p1["seed"])
        clusters = adata.obs["overclust"].astype(str)
        n_overclust = clusters.nunique()
        winner = per_cell.groupby(clusters.values).agg(lambda s: s.value_counts().index[0])
        majority = clusters.map(winner).to_numpy()
    else:
        majority = per_cell.to_numpy()

    out = pd.DataFrame(index=adata.obs_names)
    out["author_label"] = adata.obs[cols["cell_type"]].astype(str).values
    out["celltypist_label"] = per_cell.to_numpy()
    out["celltypist_majority"] = majority
    out["conf_score"] = np.concatenate(confs)
    out["author_lineage"] = out["author_label"].map(AUTHOR_LINEAGE).fillna("other")
    out["celltypist_lineage"] = out["celltypist_majority"].map(celltypist_lineage)
    out["celltypist_lineage_per_cell"] = out["celltypist_label"].map(celltypist_lineage)
    paths["reports"].mkdir(parents=True, exist_ok=True)
    out.to_csv(paths["celltypist_csv"])

    # ---- agreement table -----------------------------------------------------------
    rows = []
    for lab, g in out.groupby("author_label"):
        top = g["celltypist_majority"].value_counts(normalize=True).head(3)
        rows.append({
            "author label": lab,
            "cells": len(g),
            "agreement % (majority)": round(100 * (g["celltypist_lineage"] == g["author_lineage"]).mean(), 1),
            "agreement % (per cell)": round(100 * (g["celltypist_lineage_per_cell"] == g["author_lineage"]).mean(), 1),
            "top CellTypist labels": "; ".join(f"{k} ({v:.0%})" for k, v in top.items()),
        })
    agree = pd.DataFrame(rows).sort_values("cells", ascending=False)
    overall = 100 * (out["celltypist_lineage"] == out["author_lineage"]).mean()
    overall_cell = 100 * (out["celltypist_lineage_per_cell"] == out["author_lineage"]).mean()

    # ---- heatmap: author label x CellTypist label (row %) ------------------------------
    ct = pd.crosstab(out["author_label"], out["celltypist_majority"], normalize="index") * 100
    ct = ct.loc[agree["author label"]]
    keep = ct.columns[(ct.max(axis=0) >= 2)]
    keep = ct[keep].max(axis=0).sort_values(ascending=False).index[:20]
    ct = ct[keep]
    apply_style()
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("blue", BLUE_RAMP)
    fig, ax = plt.subplots(figsize=(0.45 * ct.shape[1] + 3.5, 0.38 * ct.shape[0] + 2.2))
    im = ax.imshow(ct.to_numpy(), cmap=cmap, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(ct.shape[1]), ct.columns, rotation=55, ha="right", fontsize=7.5)
    ax.set_yticks(range(ct.shape[0]), ct.index, fontsize=8)
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = ct.iat[i, j]
            if v >= 10:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                        color=SURFACE if v >= 55 else INK)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("% of the author label's cells", fontsize=8)
    cb.outline.set_visible(False)
    ax.set_xlabel("CellTypist label (majority vote per cluster)")
    ax.set_ylabel("Author label")
    ax.set_title("Where each author label lands in CellTypist (row %, labels covering >= 2% shown)", loc="left")
    paths["figures"].mkdir(parents=True, exist_ok=True)
    fig_path = paths["figures"] / "05_celltypist_vs_author.png"
    fig.savefig(fig_path)
    plt.close(fig)

    L = ["# Phase 1: cell-type label check with CellTypist\n"]
    vote = (f"majority voting over {n_overclust} fine Leiden clusters on the scVI graph"
            if n_overclust else "no majority voting")
    L.append(f"- Model: `{model_name}`; {vote}")
    L.append(f"- Cells: {len(out):,}; distinct CellTypist labels: {out['celltypist_majority'].nunique()}")
    L.append(f"- **Lineage agreement with the author labels: {overall:.1f}% after majority voting, "
             f"{overall_cell:.1f}% per cell**\n")
    L.append(agree.to_markdown(index=False))
    L.append(
        "\nLow agreement for a label means either the author label is coarse or ambiguous "
        "(e.g. `lymphocyte`, which has no single lineage), or the two methods genuinely "
        "disagree; the top-labels column shows which. Fine-grained CellTypist subsets "
        "(e.g. regulatory vs effector-memory helper T cells) can be used for the "
        "cell-type proportion analysis in phase 2."
    )
    L.append(f"\nFigure: `{fig_path.name}`; per-cell predictions: `{paths['celltypist_csv'].name}`")
    write_report(paths["reports"] / "05_annotation.md", L)


if __name__ == "__main__":
    main()
