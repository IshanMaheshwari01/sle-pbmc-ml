#!/usr/bin/env python
"""Phase 1, step 2: integrate across processing batches and check what it costs.

Three embeddings of the same cells are compared:
    X_pca      no batch correction (the baseline)
    X_harmony  Harmony on the PCA
    X_scVI     scVI, a variational autoencoder trained on raw counts

In this dataset batch is tangled up with disease (batch 1 is all healthy,
batch 2 almost all lupus), so batch correction can remove real disease
biology along with technical noise. Standard integration scores reward
mixing batches, which here partly means mixing lupus with healthy cells. So
this script also measures, in batch 4 (the batch with both groups), whether
cells still sit near other cells of the same disease status after each
method. Neighbours from the same donor are excluded, otherwise donor identity
alone would make neighbours share a diagnosis.

The integrated embeddings are for visualisation, clustering and annotation
only. Differential expression and the classifier (phases 2-3) use pseudobulk
raw counts with batch as a covariate, never corrected values.

Also computes a six-gene type I interferon score (Rice et al. 2013) per donor
in classical monocytes, a known lupus signal, as a sanity check that does not
depend on integration at all.

Outputs
    data/processed/<subsample>_integrated.h5ad   embeddings, UMAPs, clusters, IFN score
    results/models/<subsample>/scvi/             trained scVI model
    results/reports/04_integration.md
    results/reports/04_scib_metrics.csv
    results/figures/04_umap_integration.png, 04_umap_cell_types.png, 04_ifn_score.png

Usage:  python scripts/04_integrate.py [--config ...] [--dev] [--max-epochs N] [--skip-scib]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.phase1 import (  # noqa: E402
    gene_symbols,
    load_phase1,
    phase1_paths,
    stratified_subset,
    write_report,
)
from sle_pbmc.plotting import (  # noqa: E402
    DISEASE_COLORS,
    DISEASE_SHORT,
    INK,
    apply_style,
    plt,
    umap_highlight,
    add_disease_legend,
    umap_two_groups,
)

EMBEDDINGS = {"Unintegrated (PCA)": "X_pca", "Harmony": "X_harmony", "scVI": "X_scVI"}


# --------------------------------------------------------------------------
def run_scvi(hv: ad.AnnData, batch: str, p: dict, seed: int, model_dir: Path, max_epochs: int | None):
    import scvi
    import torch

    scvi.settings.seed = seed
    gpu = torch.cuda.is_available()
    print(f"scVI on {'GPU: ' + torch.cuda.get_device_name(0) if gpu else 'CPU (slow; see README for GPU)'}")
    scvi.model.SCVI.setup_anndata(hv, layer="counts", batch_key=batch)
    model = scvi.model.SCVI(hv, n_latent=p["n_latent"], n_layers=p["n_layers"], n_hidden=p["n_hidden"])
    t0 = time.time()
    model.train(
        max_epochs=max_epochs or p["max_epochs"],
        batch_size=p["batch_size"],
        early_stopping=True,
        accelerator="gpu" if gpu else "cpu",
        enable_progress_bar=True,
    )
    hist = model.history
    epochs = len(hist["elbo_train"]) if "elbo_train" in hist else None
    model_dir.mkdir(parents=True, exist_ok=True)
    model.save(model_dir, overwrite=True)
    return model.get_latent_representation(), {
        "device": torch.cuda.get_device_name(0) if gpu else "CPU",
        "epochs": epochs,
        "minutes": round((time.time() - t0) / 60, 1),
    }


def disease_neighbour_enrichment(
    adata: ad.AnnData, rep: str, cols: dict, batch: str, batch_value: str, k: int, min_cells: int
) -> pd.DataFrame:
    """Within one batch and one cell type: how often are a cell's nearest
    neighbours (from OTHER donors) of the same disease status, relative to chance?
    1.0 = no disease structure; higher = cells cluster by disease."""
    rows = []
    in_batch = adata.obs[batch].astype(str).values == batch_value
    for ct in adata.obs[cols["cell_type"]].astype(str).unique():
        m = in_batch & (adata.obs[cols["cell_type"]].astype(str).values == ct)
        lab = adata.obs.loc[m, cols["disease"]].astype(str).values
        donors = adata.obs.loc[m, cols["donor"]].astype(str).values
        counts = pd.Series(lab).value_counts()
        if len(counts) < 2 or counts.min() < min_cells:
            continue
        X = adata.obsm[rep][m]
        n_query = min(len(X), k * 8 + 1)
        _, idx = NearestNeighbors(n_neighbors=n_query).fit(X).kneighbors(X)
        # chance level for each cell: the share of its label among cells from OTHER donors
        tab = pd.crosstab(donors, lab)
        lab_tot = tab.sum(axis=0)
        cells_tot = int(tab.to_numpy().sum())
        same, total, exp_sum = 0, 0, 0.0
        for i in range(len(X)):
            nb = idx[i, 1:]
            nb = nb[donors[nb] != donors[i]][:k]
            same += int((lab[nb] == lab[i]).sum())
            total += len(nb)
            own = tab.loc[donors[i]]
            exp_sum += len(nb) * (lab_tot[lab[i]] - own[lab[i]]) / (cells_tot - own.sum())
        expected = exp_sum / total if total else np.nan
        observed = same / total if total else np.nan
        rows.append({"cell_type": ct, "cells": int(m.sum()), "observed": observed, "expected": expected,
                     "enrichment": observed / expected})
    return pd.DataFrame(rows)


def pick_balanced_batch(obs: pd.DataFrame, batch: str, disease: str, donor: str) -> str:
    t = obs.groupby([obs[batch].astype(str), obs[disease].astype(str)], observed=True)[donor].nunique().unstack(fill_value=0)
    return str(t.min(axis=1).idxmax())


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--phase1-config", default=None)
    ap.add_argument("--dev", action="store_true")
    ap.add_argument("--max-epochs", type=int, default=None, help="override scVI epochs (for quick tests)")
    ap.add_argument("--skip-scib", action="store_true", help="skip the scib-metrics benchmark")
    args = ap.parse_args()

    cfg, p1 = load_config(args.config), load_phase1(args.phase1_config)
    paths = phase1_paths(cfg, dev=args.dev)
    cols, batch, seed = cfg["columns"], p1["batch_key"], p1["seed"]
    if not paths["qc"].exists():
        sys.exit(f"{paths['qc']} not found. Run `make qc` first.")
    sc.settings.verbosity = 1
    t_start = time.time()

    adata = ad.read_h5ad(paths["qc"])
    adata.obs[batch] = adata.obs[batch].astype(str).astype("category")
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sym = gene_symbols(adata.var, cols["gene_symbol"])

    # ---- normalise + HVGs + PCA ------------------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_hvg = min(p1["hvg"]["n_top_genes"], adata.n_vars - 1)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_hvg, layer="counts", batch_key=batch)
    hv = adata[:, adata.var["highly_variable"].values].copy()
    n_comps = min(p1["pca"]["n_comps"], hv.n_vars - 1, hv.n_obs - 1)
    sc.pp.pca(hv, n_comps=n_comps, random_state=seed)

    # ---- Harmony ----------------------------------------------------------------
    import harmonypy

    try:
        ho = harmonypy.run_harmony(hv.obsm["X_pca"], hv.obs, [batch], random_state=seed, verbose=False)
    except TypeError:  # older/newer harmonypy without these keywords
        ho = harmonypy.run_harmony(hv.obsm["X_pca"], hv.obs, [batch])
    Z = np.asarray(ho.Z_corr)
    hv.obsm["X_harmony"] = Z.T if Z.shape[0] != hv.n_obs else Z

    # ---- scVI ---------------------------------------------------------------------
    hv.obsm["X_scVI"], scvi_info = run_scvi(hv, batch, p1["scvi"], seed, paths["models"] / "scvi", args.max_epochs)

    for key in EMBEDDINGS.values():
        adata.obsm[key] = hv.obsm[key]

    # ---- neighbours, UMAPs, clusters -----------------------------------------
    for name, key in EMBEDDINGS.items():
        tag = key.replace("X_", "")
        sc.pp.neighbors(adata, use_rep=key, n_neighbors=p1["neighbors"]["n_neighbors"], key_added=tag, random_state=seed)
        sc.tl.umap(adata, neighbors_key=tag, random_state=seed)
        adata.obsm[f"X_umap_{tag}"] = adata.obsm.pop("X_umap")
    sc.tl.leiden(adata, neighbors_key="scVI", resolution=1.0, key_added="leiden_scvi", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)

    # ---- scib-metrics benchmark ------------------------------------------------------
    scib_table, scib_note = None, ""
    if not args.skip_scib:
        try:
            from scib_metrics.benchmark import Benchmarker

            sub_idx = stratified_subset(adata.obs, [batch, cols["cell_type"]], p1["benchmark"]["n_cells"], seed)
            sub = ad.AnnData(obs=adata.obs.iloc[sub_idx][[batch, cols["cell_type"]]].copy())
            for key in EMBEDDINGS.values():
                sub.obsm[key] = np.asarray(adata.obsm[key][sub_idx], dtype=np.float32)
            sub.obs[cols["cell_type"]] = sub.obs[cols["cell_type"]].astype(str).astype("category")
            bm = Benchmarker(sub, batch_key=batch, label_key=cols["cell_type"],
                             embedding_obsm_keys=list(EMBEDDINGS.values()),
                             pre_integrated_embedding_obsm_key="X_pca", n_jobs=-1)
            bm.benchmark()
            res = bm.get_results(min_max_scale=False)
            metric_type = res.loc["Metric Type"] if "Metric Type" in res.index else None
            res = res.drop(index="Metric Type", errors="ignore").apply(pd.to_numeric, errors="coerce")
            res.index = [{v: k for k, v in EMBEDDINGS.items()}.get(i, i) for i in res.index]
            scib_table = res.T  # metrics as rows, embeddings as columns
            if metric_type is not None:
                scib_table.insert(0, "type", metric_type.reindex(scib_table.index).values)
            paths["reports"].mkdir(parents=True, exist_ok=True)
            scib_table.to_csv(paths["reports"] / "04_scib_metrics.csv")
            scib_note = f"computed on {len(sub_idx):,} cells, stratified by batch and cell type"
        except Exception as e:  # keep the rest of the analysis if scib fails
            scib_note = f"skipped: {type(e).__name__}: {e}"
            print(f"[warn] scib-metrics {scib_note}")
    else:
        scib_note = "skipped (--skip-scib)"

    # ---- disease signal retention ------------------------------------------------------
    ds = p1["disease_signal"]
    bval = str(ds["batch"])
    if bval == "auto" or bval not in set(adata.obs[batch].astype(str)):
        bval = pick_balanced_batch(adata.obs, batch, cols["disease"], cols["donor"])
    enr = []
    for name, key in EMBEDDINGS.items():
        t = disease_neighbour_enrichment(adata, key, cols, batch, bval, ds["k"], ds["min_cells_per_group"])
        t.insert(0, "embedding", name)
        enr.append(t)
    enr = pd.concat(enr, ignore_index=True) if enr else pd.DataFrame()

    # ---- interferon score ------------------------------------------------------------------
    ifn = p1["ifn"]
    sym_to_id = pd.Series(adata.var_names, index=sym)
    sym_to_id = sym_to_id[~sym_to_id.index.duplicated()]
    present = [g for g in ifn["genes"] if g in sym_to_id.index]
    missing = [g for g in ifn["genes"] if g not in sym_to_id.index]
    ifn_donor = pd.DataFrame()
    if present:
        sc.tl.score_genes(adata, gene_list=list(sym_to_id[present]), score_name="ifn_score",
                          ctrl_size=min(50, adata.n_vars // 4), random_state=seed)
        mono = adata.obs[cols["cell_type"]].astype(str) == ifn["cell_type"]
        o = adata.obs.loc[mono]
        per_sample = o.groupby(["sample_id", batch], observed=True).agg(
            ifn=("ifn_score", "median"), cells=("ifn_score", "size"),
            donor=(cols["donor"], "first"), disease=(cols["disease"], "first"),
        ).reset_index()
        per_sample = per_sample[per_sample["cells"] >= 10]
        # one value per donor per batch (a donor's flare + treated samples are averaged)
        ifn_donor = per_sample.groupby(["donor", batch], observed=True).agg(
            ifn=("ifn", "mean"), disease=("disease", "first")).reset_index()

    # ---- figures ---------------------------------------------------------------------------
    apply_style()
    paths["figures"].mkdir(parents=True, exist_ok=True)
    batches = sorted(adata.obs[batch].astype(str).unique())
    fig, axes = plt.subplots(len(EMBEDDINGS), 1 + len(batches), figsize=(2.3 * (1 + len(batches)), 2.4 * len(EMBEDDINGS)))
    axes = np.atleast_2d(axes)
    dis = adata.obs[cols["disease"]].astype(str).values
    bvals = adata.obs[batch].astype(str).values
    for r, (name, key) in enumerate(EMBEDDINGS.items()):
        xy = adata.obsm[f"X_umap_{key.replace('X_', '')}"]
        umap_two_groups(axes[r, 0], xy, dis, "Disease", seed=seed)
        for c, b in enumerate(batches, start=1):
            umap_highlight(axes[r, c], xy, bvals == b, f"Batch {b}")
        axes[r, 0].text(-0.08, 0.5, name, transform=axes[r, 0].transAxes, rotation=90,
                        va="center", ha="right", fontsize=10, color=INK)
    fig.suptitle("Same cells, three embeddings: does each batch overlap the others, and does disease still separate?",
                 fontsize=10, x=0.01, ha="left", y=1.0)
    add_disease_legend(fig)
    f1 = paths["figures"] / "04_umap_integration.png"
    fig.savefig(f1)
    plt.close(fig)

    ct_plot = adata.obs[cols["cell_type"]].astype(str).replace(p1.get("plot_cell_type_map", {}))
    groups = ct_plot.value_counts().index.tolist()
    if "Other" in groups:
        groups = [g for g in groups if g != "Other"] + ["Other"]
    ncol = 4
    nrow = int(np.ceil(len(groups) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.4 * ncol, 2.4 * nrow))
    axes = np.atleast_1d(axes).ravel()
    xy = adata.obsm["X_umap_scVI"]
    for ax, g in zip(axes, groups):
        umap_highlight(ax, xy, ct_plot.values == g, g)
    for ax in axes[len(groups):]:
        ax.axis("off")
    fig.suptitle("Author cell-type labels on the scVI UMAP", fontsize=10, x=0.01, ha="left")
    f2 = paths["figures"] / "04_umap_cell_types.png"
    fig.savefig(f2)
    plt.close(fig)

    f3 = None
    if len(ifn_donor):
        bs = sorted(ifn_donor[batch].astype(str).unique())
        fig, ax = plt.subplots(figsize=(1.6 * len(bs) + 1.5, 3.4))
        rng = np.random.default_rng(seed)
        ticks = []
        for i, b in enumerate(bs):
            ns = []
            for j, (dname, col) in enumerate(DISEASE_COLORS.items()):
                v = ifn_donor.loc[(ifn_donor[batch].astype(str) == b) & (ifn_donor["disease"] == dname), "ifn"].to_numpy()
                ns.append(f"{DISEASE_SHORT[dname]} {len(v)}")
                if not len(v):
                    continue
                x = i + (j - 0.5) * 0.36
                ax.scatter(x + rng.uniform(-0.09, 0.09, len(v)), v, s=12, c=col, alpha=0.75, linewidths=0)
                ax.plot([x - 0.13, x + 0.13], [np.median(v)] * 2, color=INK, lw=1.6)
            ticks.append(f"Batch {b}\n" + " · ".join(ns))
        ax.set_xticks(range(len(bs)), ticks, fontsize=7.5)
        ax.set_ylabel(f"IFN score in {ifn['cell_type']}s\n(median per donor)")
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        for dname, col in DISEASE_COLORS.items():
            ax.scatter([], [], s=20, c=col, label=DISEASE_SHORT[dname])
        ax.legend(loc="upper left", fontsize=8)
        ax.set_title("Type I interferon signature by batch (bar = median)", loc="left")
        f3 = paths["figures"] / "04_ifn_score.png"
        fig.savefig(f3)
        plt.close(fig)

    # ---- report ---------------------------------------------------------------------------------
    L = ["# Phase 1: integration\n"]
    L.append(f"- Cells: {adata.n_obs:,}; highly variable genes: {int(adata.var['highly_variable'].sum()):,} "
             f"(seurat_v3, ranked across `{batch}`); PCs: {n_comps}")
    L.append(f"- scVI: {scvi_info['epochs']} epochs on {scvi_info['device']} in {scvi_info['minutes']} min; "
             f"latent dimensions: {p1['scvi']['n_latent']}")
    L.append(f"- Leiden clusters on the scVI graph: {adata.obs['leiden_scvi'].nunique()}")

    L.append("\n## Standard integration benchmark (scib-metrics)\n")
    L.append(f"_{scib_note}_\n")
    if scib_table is not None:
        L.append(scib_table.round(3).to_markdown())
        L.append("\nAll metrics are 0-1, higher is better (unscaled).")
        L.append(
            "\n**Read the batch-correction scores with care.** They reward cells from different "
            "batches mixing together. Here batch 1 is all healthy and batch 2 almost all lupus, so "
            "part of what a high batch score rewards is lupus and healthy cells overlapping, i.e. "
            "lost biology. The next table checks that directly."
        )

    L.append(f"\n## Is disease signal still there? (batch {bval} only)\n")
    L.append(
        "For each cell, its nearest neighbours **from other donors** are taken in each embedding, "
        "within one cell type and within the one batch that holds both groups. *Enrichment* is how much "
        "more often those neighbours share the cell's disease status than random cells would "
        "(1.0 = none). If integration lowers enrichment compared with the unintegrated PCA, it removed "
        "disease differences along with batch."
    )
    if len(enr):
        piv = enr.pivot(index="cell_type", columns="embedding", values="enrichment")[list(EMBEDDINGS)]
        piv.insert(0, "cells", enr.groupby("cell_type")["cells"].first())
        L.append("\n" + piv.round(3).sort_values("cells", ascending=False).to_markdown())
    else:
        L.append("\n_No cell type in that batch had enough cells from both groups._")

    L.append(f"\n## Interferon score in {ifn['cell_type']}s\n")
    if missing:
        L.append(f"Genes not found and skipped: {', '.join(missing)}\n")
    if len(ifn_donor):
        L.append(
            "Six-gene type I IFN score (Rice et al. 2013) from log-normalised expression, "
            "median per sample, averaged per donor. AUROC = how well the score alone separates "
            "lupus from healthy donors (0.5 = chance). This uses no batch correction, so a high "
            "AUROC *within* a batch cannot come from batch."
        )
        rows = []
        for b, g in list(ifn_donor.groupby(batch, observed=True)) + [("all", ifn_donor)]:
            sle = g.loc[g["disease"] == "systemic lupus erythematosus", "ifn"]
            hc = g.loc[g["disease"] == "normal", "ifn"]
            row = {"batch": str(b), "SLE donors": len(sle), "healthy donors": len(hc),
                   "median SLE": sle.median() if len(sle) else np.nan,
                   "median healthy": hc.median() if len(hc) else np.nan}
            if len(sle) >= 3 and len(hc) >= 3:
                u, pval = mannwhitneyu(sle, hc, alternative="two-sided")
                row["AUROC"] = u / (len(sle) * len(hc))
                row["p (Mann-Whitney)"] = pval
            rows.append(row)
        L.append("\n" + pd.DataFrame(rows).round(4).to_markdown(index=False))
        L.append("\nThe `all` row pools batches and so mixes the disease effect with the batch effect; "
                 "the within-batch rows are the honest ones.")
    L.append("\n## Figures\n")
    for f in (f1, f2, f3):
        if f is not None:
            L.append(f"- `{f.name}`")
    L.append(f"\nTotal time: {(time.time() - t_start) / 60:.1f} min")
    write_report(paths["reports"] / "04_integration.md", L)

    adata.uns["phase1"] = {"disease_signal_batch": bval, "scvi": scvi_info}
    adata.write_h5ad(paths["integrated"], compression="gzip")
    print(f"Wrote {paths['integrated']}")


if __name__ == "__main__":
    main()
