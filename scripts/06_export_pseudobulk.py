#!/usr/bin/env python
"""Phase 2, step 1: hand the pseudobulk counts to R and define the cohort.

Writes plain files that base R reads without extra packages:
    data/processed/phase2/counts.mtx.gz   genes x profiles, raw integer counts
    data/processed/phase2/profiles.csv    one row per profile, with cohort flags
    data/processed/phase2/genes.csv       Ensembl ID + symbol
    data/processed/phase2/manifest.yaml   paths and phase-2 settings for the R scripts

A profile is one sample x processing batch x cell type (see phase 0). The
cohort rules (config/phase2.yaml) are applied here, once, so both R analyses
use exactly the same donors:
    * women only (15 of 17 men have lupus)
    * ancestry groups with enough donors
    * one biological sample per donor (a flare and a post-treatment sample of
      the same patient are not independent; the first state in
      `state_priority` is kept)

Usage:  python scripts/06_export_pseudobulk.py [--config config/config.yaml]
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.phase1 import write_report  # noqa: E402

PHASE2_CONFIG = Path(__file__).resolve().parents[1] / "config" / "phase2.yaml"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--phase2-config", default=None)
    ap.add_argument("--out-dir", default=None, help="default: <processed>/phase2")
    args = ap.parse_args()

    cfg = load_config(args.config)
    p2_path = Path(args.phase2_config) if args.phase2_config else PHASE2_CONFIG
    p2 = yaml.safe_load(p2_path.read_text())
    cols = cfg["columns"]
    pb_path = Path(cfg["pseudobulk"]["out_h5ad"])
    if not pb_path.exists():
        sys.exit(f"{pb_path} not found. Run `make subsample` first.")
    out = Path(args.out_dir) if args.out_dir else pb_path.parent / "phase2"
    out.mkdir(parents=True, exist_ok=True)

    pb = ad.read_h5ad(pb_path)
    obs = pb.obs.copy()
    for c in ("donor", "disease", "disease_state", "batch", "ancestry", "sex", "cell_type"):
        if c not in obs.columns:
            sys.exit(f"pseudobulk obs is missing `{c}`; re-run `make subsample` with the current code.")
        obs[c] = obs[c].astype(str)
    if (obs["batch"] == "mixed").any():
        sys.exit("Some profiles span several batches. Set pseudobulk.split_by: [\"Processing_Cohort\"] "
                 "in config.yaml and re-run `make subsample`.")

    # ---- cohort -------------------------------------------------------------------
    c = p2["cohort"]
    keep = obs["sex"].isin(c["sexes"])
    donors = obs.loc[keep].groupby("donor")["ancestry"].first()
    anc_counts = donors.value_counts()
    ok_anc = anc_counts[anc_counts >= c["min_ancestry_donors"]].index
    keep &= obs["ancestry"].isin(ok_anc)

    # one sample per donor: the first disease_state in the priority list
    rank = {s: i for i, s in enumerate(c["state_priority"])}
    samples = obs.loc[keep, ["donor", "sample_id", "disease_state"]].drop_duplicates("sample_id")
    samples["rank"] = samples["disease_state"].map(rank).fillna(len(rank))
    chosen = samples.sort_values(["donor", "rank"]).drop_duplicates("donor")["sample_id"]
    obs["in_cohort"] = keep & obs["sample_id"].isin(set(chosen))
    obs["disease_group"] = obs["disease"].map(
        lambda d: "healthy" if d == "normal" else ("SLE" if "lupus" in d else d))

    # ---- write ------------------------------------------------------------------------
    X = pb.X if sp.issparse(pb.X) else sp.csr_matrix(pb.X)
    with gzip.open(out / "counts.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, X.T.tocoo(), field="integer")  # genes x profiles, as R expects
    obs.insert(0, "profile", obs.index.astype(str))
    obs.to_csv(out / "profiles.csv", index=False)
    genes = pd.DataFrame({"gene_id": pb.var_names.astype(str)})
    sym_col = cols["gene_symbol"]
    genes["symbol"] = pb.var[sym_col].astype(str).values if sym_col in pb.var else genes["gene_id"]
    genes.to_csv(out / "genes.csv", index=False)

    reports = Path(cfg["reports_dir"])
    parts = list(reports.parts)
    if "reports" in parts:
        i = len(parts) - 1 - parts[::-1].index("reports")
        figures = Path(*parts[:i], "figures", *parts[i + 1:])
    else:
        figures = reports / "figures"
    manifest = {
        "counts": str((out / "counts.mtx.gz").resolve()),
        "profiles": str((out / "profiles.csv").resolve()),
        "genes": str((out / "genes.csv").resolve()),
        "reports_dir": str(reports.resolve()),
        "figures_dir": str(figures.resolve()),
        "phase2": p2,
    }
    (out / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))

    # ---- cohort report ------------------------------------------------------------------
    co = obs[obs["in_cohort"]]
    d = co.groupby("donor").agg(disease=("disease_group", "first"), ancestry=("ancestry", "first"),
                                state=("disease_state", "first"))
    L = ["# Phase 2 cohort\n"]
    L.append(f"- Profiles exported: {len(obs):,} (sample x batch x cell type), genes: {pb.n_vars:,}")
    L.append(f"- Rules: sex in {c['sexes']}; ancestry groups with >= {c['min_ancestry_donors']} donors "
             f"(kept: {', '.join(map(str, ok_anc))}); one sample per donor, priority {c['state_priority']}")
    dropped_anc = [f"{a} ({n})" for a, n in anc_counts.items() if a not in set(ok_anc)]
    if dropped_anc:
        L.append(f"- Ancestry groups left out: {', '.join(dropped_anc)}")
    L.append(f"- **Cohort: {len(d)} donors ({(d['disease'] == 'SLE').sum()} SLE, "
             f"{(d['disease'] == 'healthy').sum()} healthy), {co['profile_id'].nunique()} sample x batch "
             f"units, {len(co):,} profiles**\n")
    L.append("## Donors by ancestry\n")
    L.append(pd.crosstab(d["ancestry"], d["disease"], margins=True).to_markdown())
    L.append("\n## SLE samples by disease state\n")
    L.append(d.loc[d["disease"] == "SLE", "state"].value_counts().rename("donors").to_frame().to_markdown())
    units = co.drop_duplicates("profile_id")
    L.append("\n## Sample x batch units by batch\n")
    L.append(pd.crosstab(units["batch"], units["disease_group"], margins=True).to_markdown())
    L.append("\nA sample measured in two batches gives two units; the models treat them as correlated "
             "measurements of the same donor (limma duplicateCorrelation).")
    write_report(reports / "06_cohort.md", L)
    print(f"Wrote phase-2 inputs to {out}")


if __name__ == "__main__":
    main()
