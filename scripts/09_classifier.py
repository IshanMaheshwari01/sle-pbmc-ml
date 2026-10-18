#!/usr/bin/env python
"""Phase 3: can a blood sample tell lupus from healthy, and does it still work
when batch and ancestry cannot help?

Features (per sample x batch unit, from the phase-2 export of all 1.26M cells):
    * cell-type composition (centred log-ratio of cell counts)
    * pseudobulk log-CPM expression for seven major cell types
    * a six-gene interferon score per cell type (for the IFN-only baseline)

Models, all evaluated with the same donor-grouped, repeated 5-fold CV:
    baselines   batch + ancestry only | cell-type proportions only | IFN score only
    full        L1 logistic regression (penalty tuned in an inner grouped CV),
                random forest, gradient boosting
Gene selection, imputation and scaling are refit inside every training fold.

Stress tests, the point of the exercise:
    * AUROC pooled over all donors vs WITHIN each batch that has both groups
      (batch cannot explain a within-batch result) and within each ancestry
    * the covariate-only baseline: how much a model gets from batch alone
    * can the same features tell two batches apart among healthy donors?
    * a label-permutation null for the main model
    * does the model find lupus patients WITHOUT a high interferon signature?

Outputs: results/reports/09_classifier.md (+ csv tables), results/figures/09_*.png
Usage:   python scripts/09_classifier.py [--config ...] [--phase3-config ...] [--quick]
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import brier_score_loss, roc_auc_score, roc_curve

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.phase1 import write_report  # noqa: E402
from sle_pbmc.phase3 import (  # noqa: E402
    Dataset,
    bootstrap_auc,
    build_dataset,
    donor_level,
    ifn_status,
    model_specs,
    read_phase2,
    run_cv,
)
from sle_pbmc.plotting import INK, INK_2, MUTED, SLOT, apply_style, plt  # noqa: E402

PHASE3_CONFIG = Path(__file__).resolve().parents[1] / "config" / "phase3.yaml"
GREY = "#a9a8a2"


def figures_dir(reports: Path) -> Path:
    parts = list(reports.parts)
    if "reports" in parts:
        i = len(parts) - 1 - parts[::-1].index("reports")
        return Path(*parts[:i], "figures", *parts[i + 1:])
    return reports / "figures"


def subset(ds: Dataset, mask: np.ndarray) -> Dataset:
    return Dataset(X=ds.X[mask], meta=ds.meta[mask], y=ds.y[mask], blocks=ds.blocks)


def fmt_ci(v, lo, hi):
    return "" if np.isnan(v) else f"{v:.2f} [{lo:.2f}-{hi:.2f}]"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--phase3-config", default=None)
    ap.add_argument("--data-dir", default=None, help="phase-2 export folder (default: <processed>/phase2)")
    ap.add_argument("--quick", action="store_true", help="1 repeat, 10 permutations: a fast trial run")
    args = ap.parse_args()

    t_start = time.time()
    cfg = load_config(args.config)
    p3 = yaml.safe_load(Path(args.phase3_config or PHASE3_CONFIG).read_text())
    if args.quick:
        p3["cv"]["repeats"] = 1
        p3["evaluation"]["permutations"] = min(10, p3["evaluation"]["permutations"])
    seed, ev = p3["seed"], p3["evaluation"]
    data_dir = Path(args.data_dir) if args.data_dir else Path(cfg["pseudobulk"]["out_h5ad"]).parent / "phase2"
    if not (data_dir / "counts.mtx.gz").exists():
        sys.exit(f"{data_dir} has no phase-2 export. Run `make export` first.")
    reports = Path(cfg["reports_dir"])
    figs = figures_dir(reports)
    reports.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    print("Building features ...")
    counts, profiles, genes = read_phase2(data_dir)
    ds = build_dataset(counts, profiles, genes, p3)
    meta = ds.meta
    n_donors = meta["donor"].nunique()
    print(f"{len(ds.y)} units from {n_donors} donors; {ds.X.shape[1]:,} feature columns")

    # ---- cross-validated predictions for every model --------------------------------------
    specs = model_specs(ds)
    main_name = ev["importance_model"]
    results = {}
    for spec in specs:
        t0 = time.time()
        results[spec.name] = run_cv(spec, ds, p3, seed=seed, importance=(spec.name == main_name))
        print(f"  done in {time.time() - t0:.0f}s")
    avg = {k: np.nanmean(r.preds, axis=0) for k, r in results.items()}

    # ---- evaluation settings ------------------------------------------------------------------
    settings = [("All donors (pooled)", None, "donor")]
    for b in sorted(meta["batch"].unique(), key=lambda s: float(s) if s.replace(".", "").isdigit() else s):
        m = (meta["batch"] == b).to_numpy()
        n_pos, n_neg = int(ds.y[m].sum()), int((1 - ds.y[m]).sum())
        if min(n_pos, n_neg) >= ev["min_per_class"]:
            settings.append((f"Batch {b} only", m, "unit"))
    for a in sorted(meta["ancestry"].unique()):
        m = (meta["ancestry"] == a).to_numpy()
        settings.append((f"{a} donors only", m, "donor"))

    rows, per_setting_counts = [], {}
    for name, mask, level in settings:
        for spec in specs:
            p = avg[spec.name]
            if level == "donor":
                d = donor_level(meta if mask is None else meta[mask], ds.y if mask is None else ds.y[mask],
                                p if mask is None else p[mask])
                yy, pp = d["y"].to_numpy(), d["p"].to_numpy()
            else:
                yy, pp = ds.y[mask], p[mask]
            n_pos, n_neg = int(yy.sum()), int(len(yy) - yy.sum())
            per_setting_counts[name] = (n_pos, n_neg)
            if min(n_pos, n_neg) < ev["min_per_class"]:
                continue
            if level == "unit" and spec.name.startswith("Covariates"):
                # Within a batch the covariate model has (almost) nothing to go on, and
                # refitting it fold by fold makes pooled held-out predictions anti-correlated
                # with the labels, so its AUROC there is an artefact, not information.
                continue
            auc = roc_auc_score(yy, pp)
            lo, hi = bootstrap_auc(yy, pp, ev["bootstrap"], seed)
            rows.append({"setting": name, "model": spec.name, "baseline": spec.baseline,
                         "SLE": n_pos, "healthy": n_neg, "AUROC": auc, "ci_low": lo, "ci_high": hi,
                         "brier": brier_score_loss(yy, pp)})
    perf = pd.DataFrame(rows)
    perf.to_csv(reports / "09_performance.csv", index=False)

    # spread across CV repeats (pooled, donor level)
    rep_sd = {}
    for spec in specs:
        vals = []
        for r_ in results[spec.name].preds:
            d = donor_level(meta, ds.y, r_)
            vals.append(roc_auc_score(d["y"], d["p"]))
        rep_sd[spec.name] = (np.mean(vals), np.std(vals))

    # ---- batch detectability among healthy donors -----------------------------------------------
    print("Checking how well the features identify batch among healthy donors ...")
    healthy = meta["disease_group"] == "healthy"
    hb = meta.loc[healthy, "batch"].value_counts()
    batch_check = None
    if len(hb) >= 2 and hb.iloc[1] >= ev["min_per_class"]:
        b1, b2 = hb.index[:2]
        m = (healthy & meta["batch"].isin([b1, b2])).to_numpy()
        sub = subset(ds, m)
        yb = (sub.meta["batch"] == b1).astype(int).to_numpy()
        main_spec = next(s for s in specs if s.name == main_name)
        rb = run_cv(main_spec, sub, p3, seed=seed, repeats=1, y=yb, verbose=False)  # penalty tuned for this task
        batch_check = (b1, b2, int(yb.sum()), int(len(yb) - yb.sum()), roc_auc_score(yb, rb.preds[0]))

    # ---- permutation null for the main model ---------------------------------------------------------
    print(f"Label-permutation null ({ev['permutations']} permutations) ...")
    main_spec = next(s for s in specs if s.name == main_name)
    C_mode = Counter(results[main_name].chosen_C).most_common(1)[0][0] if results[main_name].chosen_C else 1.0
    best_batch = None
    batch_settings = [s for s in settings if s[2] == "unit"]
    if batch_settings:
        best_batch = max(batch_settings, key=lambda s: min(ds.y[s[1]].sum(), (1 - ds.y[s[1]]).sum()))
    rng = np.random.default_rng(seed + 7)
    donors = meta["donor"].to_numpy()
    donor_y = pd.Series(ds.y, index=donors).groupby(level=0).first()
    null_pooled, null_batch = [], []
    for i in range(ev["permutations"]):
        perm = pd.Series(rng.permutation(donor_y.to_numpy()), index=donor_y.index)
        yp = perm.loc[donors].to_numpy()
        rp = run_cv(main_spec, ds, p3, seed=seed + 1000 + i, repeats=1, y=yp, C_fixed=C_mode, verbose=False)
        d = donor_level(meta, yp, rp.preds[0])
        null_pooled.append(roc_auc_score(d["y"], d["p"]))
        if best_batch is not None:
            mb = best_batch[1]
            null_batch.append(roc_auc_score(yp[mb], rp.preds[0][mb]) if len(np.unique(yp[mb])) == 2 else np.nan)
    obs_pooled = perf.query("setting == 'All donors (pooled)' and model == @main_name")["AUROC"].iloc[0]
    p_pooled = (1 + sum(v >= obs_pooled for v in null_pooled)) / (1 + len(null_pooled))
    obs_batch = p_batch = np.nan
    if best_batch is not None:
        obs_batch = perf.query("setting == @best_batch[0] and model == @main_name")["AUROC"].iloc[0]
        nb = [v for v in null_batch if not np.isnan(v)]
        p_batch = (1 + sum(v >= obs_batch for v in nb)) / (1 + len(nb))

    # ---- interferon-low patients -------------------------------------------------------------------------
    status = ifn_status(ds, ev["ifn_high_z"])
    ifn_rows = []
    for spec in specs:
        p = avg[spec.name]
        for grp in ("SLE IFN-high", "SLE IFN-low"):
            for label, mask in [("pooled", np.ones(len(p), bool))] + (
                    [(best_batch[0], best_batch[1])] if best_batch is not None else []):
                sel = mask & status.isin([grp, "healthy"]).to_numpy()
                yy = (status[sel] != "healthy").astype(int).to_numpy()
                if min(yy.sum(), len(yy) - yy.sum()) < ev["min_per_class"]:
                    continue
                ifn_rows.append({"model": spec.name, "SLE group": grp, "setting": label,
                                 "SLE units": int(yy.sum()), "healthy units": int(len(yy) - yy.sum()),
                                 "AUROC vs healthy": roc_auc_score(yy, p[sel])})
    ifn_tab = pd.DataFrame(ifn_rows)
    sle_units = status[status.str.startswith("SLE")]
    sle_donor_status = (pd.DataFrame({"donor": meta["donor"], "s": status})
                        .loc[meta["disease_group"] == "SLE"].groupby("donor")["s"].first())

    # ---- importance -----------------------------------------------------------------------------------------
    imp = pd.DataFrame(results[main_name].importance)
    imp_summary = imp.agg(["mean", "std"]).T.sort_values("mean", ascending=False) if len(imp) else pd.DataFrame()
    coef_tab = pd.DataFrame()
    if results[main_name].coefs:
        C = pd.concat(results[main_name].coefs, axis=1).fillna(0.0)
        nz = (C != 0)
        coef_tab = pd.DataFrame({"selected in % of folds": 100 * nz.mean(axis=1),
                                 "mean coefficient": C.where(nz).mean(axis=1)})
        coef_tab = coef_tab[coef_tab["selected in % of folds"] > 0].sort_values(
            ["selected in % of folds", "mean coefficient"], ascending=[False, False])
        coef_tab.index = [i.replace("expr|", "").replace("prop|", "proportion: ") for i in coef_tab.index]
        coef_tab.to_csv(reports / "09_selected_features.csv")

    # ---- figures ----------------------------------------------------------------------------------------------
    apply_style()
    model_order = [s.name for s in specs]
    set_names = [s[0] for s in settings if s[0] in set(perf["setting"])]
    fig, axes = plt.subplots(1, len(set_names), figsize=(2.55 * len(set_names) + 2.6, 0.42 * len(specs) + 1.6),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, sname in zip(axes, set_names):
        sub = perf[perf["setting"] == sname].set_index("model").reindex(model_order)
        ypos = np.arange(len(model_order))[::-1]
        for yv, (mname, r) in zip(ypos, sub.iterrows()):
            if pd.isna(r["AUROC"]):
                continue
            col = GREY if r["baseline"] else SLOT[0]
            ax.plot([r["ci_low"], r["ci_high"]], [yv, yv], color=col, lw=2, solid_capstyle="round")
            ax.scatter([r["AUROC"]], [yv], s=36, color=col, zorder=3, edgecolor="#fcfcfb", linewidth=1.2)
            ax.text(r["ci_high"] + 0.015, yv, f"{r['AUROC']:.2f}", va="center", fontsize=7.5, color=INK_2)
        ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)))
        ax.set_xlim(min(0.3, float(np.nanmin(perf["ci_low"])) - 0.05), 1.08)
        ax.set_xticks([0.4, 0.6, 0.8, 1.0])
        n_pos, n_neg = per_setting_counts.get(sname, (0, 0))
        ax.set_title(f"{sname}\n{n_pos} SLE vs {n_neg} healthy", fontsize=8.5)
        ax.grid(axis="x")
        ax.set_axisbelow(True)
        ax.set_xlabel("AUROC (95% CI)")
    axes[0].set_yticks(np.arange(len(model_order))[::-1], model_order)
    fig.suptitle("Lupus vs healthy: accuracy pooled, and where batch or ancestry cannot help "
                 "(grey = baselines, dashed = chance)", fontsize=10, x=0.01, ha="left", y=1.02)
    f1 = figs / "09_auroc_by_setting.png"
    fig.savefig(f1)
    plt.close(fig)

    f2 = None
    if best_batch is not None:
        mb = best_batch[1]
        fig, ax = plt.subplots(figsize=(4.8, 4.4))
        curves = [(main_name, SLOT[0]), ("Interferon score only", SLOT[1]),
                  ("Cell-type proportions only", SLOT[2])]
        for mname, col in curves:
            fpr, tpr, _ = roc_curve(ds.y[mb], avg[mname][mb])
            auc = roc_auc_score(ds.y[mb], avg[mname][mb])
            ax.plot(fpr, tpr, color=col, lw=2, label=f"{mname} ({auc:.2f})")
        ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls=(0, (3, 3)))
        ax.set_xlabel("False positive rate (healthy called SLE)")
        ax.set_ylabel("True positive rate (SLE found)")
        n_pos, n_neg = per_setting_counts[best_batch[0]]
        ax.set_title(f"ROC within {best_batch[0].replace(' only', '')} ({n_pos} SLE, {n_neg} healthy)", loc="left")
        ax.legend(loc="lower right", fontsize=7.5)
        ax.set_aspect("equal")
        ax.grid(True)
        f2 = figs / "09_roc_within_batch.png"
        fig.savefig(f2)
        plt.close(fig)

    f3 = None
    if len(imp_summary):
        fig, ax = plt.subplots(figsize=(6.2, 0.36 * len(imp_summary) + 1.4))
        lab = [s.replace(", alpha-beta T cell", " T cell").replace("-positive", "+") for s in imp_summary.index]
        ypos = np.arange(len(lab))[::-1]
        ax.barh(ypos, imp_summary["mean"], color=SLOT[0], height=0.6)
        ax.errorbar(imp_summary["mean"], ypos, xerr=imp_summary["std"], fmt="none", ecolor=INK_2, lw=0.8)
        ax.set_yticks(ypos, lab)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_xlabel("Drop in held-out AUROC when the block is shuffled (mean ± SD over folds)")
        ax.set_title(f"What the {main_name.split(' (')[0]} model relies on", loc="left")
        ax.grid(axis="x")
        ax.set_axisbelow(True)
        f3 = figs / "09_feature_block_importance.png"
        fig.savefig(f3)
        plt.close(fig)

    order = ["healthy", "SLE IFN-low", "SLE IFN-high"]
    fig, ax = plt.subplots(figsize=(5, 3.6))
    p = avg[main_name]
    rngj = np.random.default_rng(seed)
    for i, g in enumerate(order):
        v = p[(status == g).to_numpy()]
        if not len(v):
            continue
        col = SLOT[0] if g == "healthy" else SLOT[1]
        ax.scatter(i + rngj.uniform(-0.18, 0.18, len(v)), v, s=11, color=col, alpha=0.7, linewidths=0)
        ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color=INK, lw=1.8)
    ax.set_xticks(range(3), [f"{g}\n({int((status == g).sum())} units)" for g in order])
    ax.set_ylabel("Predicted probability of SLE\n(held-out, averaged over CV repeats)")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    ax.set_title("Does the model find lupus without a high interferon signature?", loc="left")
    f4 = figs / "09_predictions_by_ifn_status.png"
    fig.savefig(f4)
    plt.close(fig)

    # ---- report ------------------------------------------------------------------------------------------------------
    L = ["# Phase 3: donor-level lupus classifier\n"]
    L.append(f"- Cohort: {n_donors} donors ({int(donor_level(meta, ds.y, ds.y.astype(float))['y'].sum())} SLE), "
             f"{len(ds.y)} sample x batch units; same cohort as phase 2")
    blk = {k: len(v) for k, v in ds.blocks.items()}
    L.append(f"- Features: {blk.get('cell-type proportions', 0)} cell-type proportions (CLR); "
             f"{sum(v for k, v in blk.items() if k.endswith('expression')):,} expressed genes across "
             f"{sum(1 for k in blk if k.endswith('expression'))} cell types, of which the "
             f"{p3['features']['top_genes_per_block']} most variable per cell type are chosen inside each training fold")
    L.append(f"- CV: {p3['cv']['outer_folds']}-fold, grouped by donor, stratified, {p3['cv']['repeats']} repeats; "
             f"L1 penalty tuned in an inner {p3['cv']['inner_folds']}-fold grouped CV "
             f"(chosen C: {dict(Counter(results[main_name].chosen_C))})")
    L.append("- AUROC: 0.5 = chance, 1.0 = perfect. Brackets: 95% bootstrap CI over donors.\n")

    L.append("## Headline: accuracy pooled vs where batch and ancestry cannot help\n")
    wide = perf.assign(cell=[fmt_ci(a, b, c) for a, b, c in zip(perf.AUROC, perf.ci_low, perf.ci_high)])
    wide = wide.pivot(index="model", columns="setting", values="cell").reindex(model_order)[set_names]
    wide = wide.fillna("n/a")
    L.append(wide.to_markdown())
    L.append("\n`n/a`: the covariate-only model is not scored within a batch, where batch is constant; "
             "refitting a model with no usable signal fold by fold pushes its pooled held-out AUROC "
             "below 0.5, an artefact rather than information.")
    L.append("\nRead across a row: the pooled number can lean on batch (batch 1 is all healthy, batch 2 "
             "mostly lupus). Within one batch, batch is constant and cannot separate anyone, and within "
             "one ancestry group the same holds for ancestry. A model that keeps its accuracy in those "
             "columns has learned biology; the covariate-only row shows what batch and ancestry give for free.")

    L.append("\n## Stability across CV repeats (pooled, donor level)\n")
    L.append(pd.DataFrame({"mean AUROC": {k: v[0] for k, v in rep_sd.items()},
                           "SD over repeats": {k: v[1] for k, v in rep_sd.items()}}).reindex(model_order)
             .round(3).to_markdown())

    if batch_check:
        b1, b2, n1, n2, auc_b = batch_check
        L.append("\n## How much batch is written into the features?\n")
        L.append(f"Among **healthy donors only**, the main model's features separate batch {b1} ({n1} units) "
                 f"from batch {b2} ({n2} units) with AUROC **{auc_b:.2f}** (same grouped CV). "
                 + ("Batch is strongly encoded in expression, so pooled accuracy alone would overstate what the "
                    "model knows about lupus." if auc_b > 0.8 else
                    "Batch leaves only a modest mark on these features."))

    L.append("\n## Permutation test (main model)\n")
    L.append(f"Disease labels were shuffled across donors {len(null_pooled)} times and the whole CV re-run "
             f"(penalty fixed at C = {C_mode}).")
    L.append(f"- Pooled: observed AUROC {obs_pooled:.2f}; null mean {np.mean(null_pooled):.2f} "
             f"(95th percentile {np.percentile(null_pooled, 95):.2f}); permutation p = {p_pooled:.3g}")
    if best_batch is not None:
        L.append(f"- {best_batch[0]}: observed {obs_batch:.2f}; null mean {np.nanmean(null_batch):.2f} "
                 f"(95th percentile {np.nanpercentile(null_batch, 95):.2f}); p = {p_batch:.3g}")
    L.append(f"- With {len(null_pooled)} permutations the smallest possible p is {1 / (1 + len(null_pooled)):.3g}.")

    L.append("\n## Lupus patients without a high interferon signature\n")
    n_hi = int((sle_donor_status == "SLE IFN-high").sum())
    n_lo = int((sle_donor_status == "SLE IFN-low").sum())
    L.append(f"An SLE donor counts as IFN-high when the classical-monocyte IFN score is more than "
             f"{ev['ifn_high_z']} SD above healthy donors of the same batch: **{n_hi} IFN-high, {n_lo} IFN-low** "
             f"SLE donors. The table shows how well each model separates each group from healthy donors.")
    if len(ifn_tab):
        t = ifn_tab.pivot_table(index="model", columns=["setting", "SLE group"], values="AUROC vs healthy")
        t = t.reindex(model_order)
        t.columns = [f"{a}: {b}" for a, b in t.columns]
        L.append("\n" + t.round(2).to_markdown())
        L.append("\nThe IFN-only baseline should be near-perfect for IFN-high patients by construction. The "
                 "IFN-low columns are the real test: anything above chance there comes from biology other than "
                 "the interferon signature.")

    if len(imp_summary):
        L.append(f"\n## What the main model relies on\n")
        L.append("Held-out AUROC drop when one block of features is shuffled (larger = more relied on):\n")
        L.append(imp_summary.round(3).to_markdown())
    if len(coef_tab):
        L.append("\nFeatures most consistently selected by the L1 model across folds "
                 f"(of {len(results[main_name].coefs)} fits; positive coefficient = higher in SLE):\n")
        L.append(coef_tab.head(25).round(3).to_markdown())

    L.append("\n## Figures\n")
    for f in (f1, f2, f3, f4):
        if f is not None:
            L.append(f"- `{f.name}`")
    L.append(f"\nTotal time: {(time.time() - t_start) / 60:.1f} min")
    write_report(reports / "09_classifier.md", L)


if __name__ == "__main__":
    main()
