#!/usr/bin/env python
"""Phase 4: are lupus and healthy samples separated by 10x pool as well as by batch?

This is a multiplexed design: each 10x library (`library_uuid`) holds cells
from several donors. If lupus and healthy samples tended to sit in different
libraries within a batch, library-level technical differences could make the
"within-batch" classifier scores look better than the biology deserves.

Two questions, answered per processing batch:
  1. Design: how mixed are the libraries? (share of libraries holding both
     groups; how well a unit's library mates predict its diagnosis)
  2. Model: the WITHIN-LIBRARY AUROC of the main classifier: the chance that
     a lupus sample scores above a healthy sample run in the SAME library.
     Only pairs that share a library are compared, so neither batch nor
     library can explain the result.

The classifier is re-run here with the phase-3 settings (fewer repeats by
default) and its held-out predictions are saved.

Outputs: results/reports/10_pool_check.md, 10_unit_predictions.csv; the
"Pool check" section of README.md is filled in (between its markers).
Usage:   python scripts/10_pool_check.py [--config ...] [--repeats 2] [--no-readme]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.h5ad_stream import read_obs  # noqa: E402
from sle_pbmc.phase1 import write_report  # noqa: E402
from sle_pbmc.phase3 import build_dataset, model_specs, read_phase2, run_cv  # noqa: E402
from sle_pbmc.sampling import make_sample_ids  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MAIN = "L1 logistic (all features)"
START, END = "<!-- POOL_CHECK_START -->", "<!-- POOL_CHECK_END -->"
CRAMER = "Cramer's V (library x diagnosis)"


def within_library_auc(units: pd.DataFrame, p: pd.Series, libs: pd.Series) -> tuple[float, int, list]:
    """Concordance over (SLE, healthy) pairs that share a library.
    Returns AUROC, number of pairs and the per-library (concordant, pairs) list."""
    per_lib = []
    lib_to_units: dict[str, list[str]] = {}
    for u, ls in libs.items():
        for lib in ls:
            lib_to_units.setdefault(lib, []).append(u)
    for lib, us in lib_to_units.items():
        g = units.loc[us, "disease_group"]
        sle = p[g[g == "SLE"].index].to_numpy()
        hc = p[g[g == "healthy"].index].to_numpy()
        if len(sle) and len(hc):
            diff = sle[:, None] - hc[None, :]
            conc = float((diff > 0).sum() + 0.5 * (diff == 0).sum())
            per_lib.append((lib, conc, diff.size))
    n = sum(x[2] for x in per_lib)
    return (sum(x[1] for x in per_lib) / n if n else np.nan), n, per_lib


def bootstrap_libs(per_lib: list, n: int, seed: int) -> tuple[float, float]:
    if len(per_lib) < 3:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    conc = np.array([x[1] for x in per_lib])
    pairs = np.array([x[2] for x in per_lib])
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(per_lib), len(per_lib))
        vals.append(conc[idx].sum() / pairs[idx].sum())
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def library_mixing(units: pd.DataFrame, libs: pd.Series) -> tuple[float, float]:
    """Two plain measures of how mixed the libraries are.
    * share of units that have at least one unit of the OTHER group in one of their libraries
      (only those units can enter a within-library comparison)
    * Cramer's V between library and diagnosis over (unit, library) memberships:
      0 = diagnosis independent of library, 1 = every library holds a single group."""
    pairs = [(u, lib, units.loc[u, "disease_group"]) for u, ls in libs.items() for lib in ls]
    df = pd.DataFrame(pairs, columns=["unit", "lib", "group"])
    groups_in_lib = df.groupby("lib")["group"].agg(set)
    partnered = df.assign(other=[len(groups_in_lib[l] - {g}) > 0 for l, g in zip(df["lib"], df["group"])])
    pct = 100 * partnered.groupby("unit")["other"].any().mean()
    tab = pd.crosstab(df["lib"], df["group"])
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return pct, np.nan
    n = tab.to_numpy().sum()
    expected = np.outer(tab.sum(1), tab.sum(0)) / n
    chi2 = ((tab.to_numpy() - expected) ** 2 / expected).sum()
    v = np.sqrt(chi2 / (n * (min(tab.shape) - 1)))
    return float(pct), float(v)


def verdict(auc: float) -> str:
    if np.isnan(auc):
        return "could not be computed (no library holds both groups)"
    if auc >= 0.85:
        return "holds: the model separates lupus from healthy samples run in the same 10x library"
    if auc >= 0.7:
        return ("drops but stays well above chance: part of the within-batch accuracy may come from "
                "library differences, but most of the signal survives")
    return ("largely disappears: library (pool) differences explain much of the within-batch "
            "accuracy, so the within-batch numbers should not be read as biology")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--phase3-config", default=None)
    ap.add_argument("--repeats", type=int, default=2, help="CV repeats for the re-run classifier")
    ap.add_argument("--no-readme", action="store_true", help="do not edit README.md")
    ap.add_argument("--readme", default=str(ROOT / "README.md"))
    ap.add_argument("--summary", default=str(ROOT / "SUMMARY.md"))
    args = ap.parse_args()

    cfg = load_config(args.config)
    p3 = yaml.safe_load(Path(args.phase3_config or ROOT / "config/phase3.yaml").read_text())
    cols = cfg["columns"]
    raw = Path(cfg["data"]["raw_h5ad"])
    data_dir = Path(cfg["pseudobulk"]["out_h5ad"]).parent / "phase2"
    reports = Path(cfg["reports_dir"])
    for p in (raw, data_dir / "profiles.csv"):
        if not p.exists():
            sys.exit(f"{p} not found (needs `make download` and `make export`).")

    # ---- which libraries hold which sample x batch units -------------------------------
    print("Reading cell metadata ...")
    obs = read_obs(raw)
    split_by = list(cfg["pseudobulk"].get("split_by") or [])
    unit_id = make_sample_ids(obs, list(cfg["sample_key"]) + split_by)
    counts, profiles, genes = read_phase2(data_dir)
    cohort = set(profiles.loc[profiles["in_cohort"], "profile_id"])
    keep = unit_id.isin(cohort).to_numpy()
    libs = (pd.DataFrame({"unit": unit_id[keep].to_numpy(), "lib": obs.loc[keep, cols["library"]].astype(str).to_numpy()})
              .groupby("unit")["lib"].agg(lambda s: sorted(set(s))))

    # ---- classifier predictions -------------------------------------------------------------
    print(f"Re-running the main classifier ({args.repeats} CV repeats) ...")
    ds = build_dataset(counts, profiles, genes, p3)
    specs = {s.name: s for s in model_specs(ds)}
    preds = {}
    for name in (MAIN, "Interferon score only", "Cell-type proportions only"):
        r = run_cv(specs[name], ds, p3, seed=p3["seed"], repeats=args.repeats, verbose=(name == MAIN))
        preds[name] = pd.Series(np.nanmean(r.preds, axis=0), index=ds.X.index)
    units = ds.meta.copy()
    libs = libs.reindex(units.index)
    if libs.isna().any():
        sys.exit(f"{int(libs.isna().sum())} cohort units have no library; check `columns.library` in config.yaml.")
    out = units[["donor", "batch", "ancestry", "disease_group"]].copy()
    for k, v in preds.items():
        out[f"p_{k}"] = v
    out["libraries"] = libs.map(lambda l: ";".join(l))
    reports.mkdir(parents=True, exist_ok=True)
    out.to_csv(reports / "10_unit_predictions.csv")

    # ---- per batch -------------------------------------------------------------------------------------
    rows = []
    groups = [(f"Batch {b}", units.index[units["batch"] == b]) for b in sorted(units["batch"].unique())]
    groups.append(("All batches", units.index))
    for name, idx in groups:
        u = units.loc[idx]
        l = libs.loc[idx]
        lib_sets: dict[str, set] = {}
        for unit, ls in l.items():
            for lib in ls:
                lib_sets.setdefault(lib, set()).add(u.loc[unit, "disease_group"])
        n_libs = len(lib_sets)
        n_mixed = sum(len(v) == 2 for v in lib_sets.values())
        row = {"setting": name, "SLE units": int((u["disease_group"] == "SLE").sum()),
               "healthy units": int((u["disease_group"] == "healthy").sum()),
               "libraries": n_libs, "libraries with both groups": n_mixed}
        row["% units with an opposite-group partner"], row[CRAMER] = library_mixing(u, l)
        for k, p in preds.items():
            auc, n_pairs, per_lib = within_library_auc(u, p.loc[idx], l)
            lo, hi = bootstrap_libs(per_lib, 1000, p3["seed"])
            short = {MAIN: "model", "Interferon score only": "IFN only", "Cell-type proportions only": "proportions"}[k]
            row[f"within-library AUROC ({short})"] = auc
            if k == MAIN:
                row["pairs compared"] = n_pairs
                row["model 95% CI"] = "" if np.isnan(lo) else f"{lo:.2f}-{hi:.2f}"
            yy = (u["disease_group"] == "SLE").astype(int).to_numpy()
            if short == "model":
                row["ordinary AUROC (model)"] = roc_auc_score(yy, p.loc[idx]) if len(set(yy)) == 2 else np.nan
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(reports / "10_pool_check.csv", index=False)

    overall = tab.set_index("setting").loc["All batches"]
    w = overall["within-library AUROC (model)"]
    L = ["# Phase 4: pool (10x library) check\n",
         f"- Cohort units: {len(units)} (sample x batch), libraries: {overall['libraries']}, "
         f"of which {overall['libraries with both groups']} hold both lupus and healthy samples",
         f"- Classifier: {MAIN}, re-run with {args.repeats} CV repeats (same features and folds as phase 3)\n",
         "## Per batch\n", tab.round(3).to_markdown(index=False), "",
         "- **% units with an opposite-group partner**: share of samples that sit in a library with at "
         "least one sample of the other group. Only these can enter a within-library comparison.",
         "- **Cramer's V (library x diagnosis)**: 0 = diagnosis unrelated to library (well mixed); "
         "1 = every library holds only lupus or only healthy samples.",
         "- **within-library AUROC**: the chance that a lupus sample scores above a healthy sample "
         "**from the same library**. Batch and library are identical within each comparison, so they "
         "cannot explain it.",
         "- **ordinary AUROC**: the usual score over all units in that setting, for comparison.\n",
         f"## Verdict\n\nWithin-library AUROC over all batches: **{w:.2f}** "
         f"({overall['pairs compared']:,} lupus-healthy pairs sharing a library; 95% CI over libraries "
         f"{overall['model 95% CI']}). The classifier's signal {verdict(w)}."]
    write_report(reports / "10_pool_check.md", L)

    if not args.no_readme:
        readme = Path(args.readme)
        if readme.exists() and START in readme.read_text():
            b = tab.set_index("setting")
            lines = ["| Setting | Libraries (with both groups) | Samples with a same-library partner | Cramér's V, library × diagnosis | Within-library AUROC, model | Within-library, IFN only | Ordinary AUROC, model |",
                     "|---|---|---|---|---|---|---|"]
            for s, r in b.iterrows():
                f = lambda v: "" if pd.isna(v) else f"{v:.2f}"  # noqa: E731
                cv = f(r[CRAMER])
                lines.append(f"| {str(s).replace('.0', '')} | {r['libraries']} ({r['libraries with both groups']}) | "
                             f"{r['% units with an opposite-group partner']:.0f}% | {cv} | "
                             f"{f(r['within-library AUROC (model)'])} | {f(r['within-library AUROC (IFN only)'])} | "
                             f"{f(r['ordinary AUROC (model)'])} |")
            block = "\n".join([START, "",
                               f"Comparing only lupus and healthy samples run in the same 10x library "
                               f"(so neither batch nor library can help), the classifier's AUROC is {w:.2f} "
                               f"(95% CI {overall['model 95% CI']}, {overall['pairs compared']:,} pairs). "
                               f"The signal {verdict(w)}.", "", *lines, "",
                               "This section is filled in by `make poolcheck`. Details are in `results/reports/10_pool_check.md`.",
                               "", END])
            txt = re.sub(re.escape(START) + r".*?" + re.escape(END), lambda _: block, readme.read_text(), flags=re.S)
            readme.write_text(txt)
            print(f"Updated the pool-check section of {readme}")
        summary = Path(args.summary)
        s0, s1 = "<!-- POOL_SENTENCE_START -->", "<!-- POOL_SENTENCE_END -->"
        if summary.exists() and s0 in summary.read_text():
            sentence = (f"Same batch, same 10x library: comparing only SLE and healthy samples run in the "
                        f"same library, AUROC is {w:.2f} (95% CI {overall['model 95% CI']}), so the signal "
                        f"{verdict(w).split(':')[0]}.")
            txt = re.sub(re.escape(s0) + r".*?" + re.escape(s1), lambda _: s0 + sentence + s1,
                         summary.read_text(), flags=re.S)
            summary.write_text(txt)
            print(f"Updated the pool-check line of {summary}")


if __name__ == "__main__":
    main()
