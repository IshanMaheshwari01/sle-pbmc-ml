#!/usr/bin/env python
"""Audit the raw h5ad before any analysis: structure, design and confounding.

Reads only the cell metadata (obs), gene table and a slice of the count matrix,
so it runs in a minute or two even on the 12 GB file. Writes
results/reports/01_data_audit.md and prints the same tables.

The point of this step is to know, before modelling, how disease status lines
up with processing batch, ancestry and sex. If a covariate is unevenly spread
across cases and controls, a classifier can learn it instead of the disease.

Usage:  python scripts/01_inspect.py [--config config/config.yaml]
"""

from __future__ import annotations

import argparse
import difflib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402
from sle_pbmc.h5ad_stream import (  # noqa: E402
    check_integer_counts,
    locate_counts,
    read_obs,
    read_var,
)
from sle_pbmc.sampling import make_sample_ids  # noqa: E402


def md_table(df: pd.DataFrame | pd.Series) -> str:
    if isinstance(df, pd.Series):
        df = df.to_frame()
    return df.to_markdown()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    path = Path(cfg["data"]["raw_h5ad"])
    cols = cfg["columns"]
    if not path.exists():
        sys.exit(f"{path} not found. Run `make download` first.")

    obs = read_obs(path)
    loc = locate_counts(path)
    var = read_var(path, loc)
    counts = check_integer_counts(path, loc)

    out: list[str] = ["# Data audit: Perez et al. 2022 lupus PBMC\n"]
    out.append(f"- File: `{path.name}`")
    out.append(f"- Cells: **{len(obs):,}**; genes in count matrix: **{loc.shape[1]:,}**")
    out.append(f"- Raw counts read from: `{loc.matrix_path}`")
    out.append(
        f"- First {counts['checked_values']:,} stored values: "
        f"{'all whole numbers' if counts['all_integer'] else 'NOT all whole numbers'} "
        f"(min {counts['min']}, max {counts['max']}, stored as {counts['dtype_on_disk']})"
    )
    if not counts["all_integer"]:
        out.append(
            "\n> **Warning:** these do not look like raw counts. DESeq2/edgeR need "
            "raw integer counts; check which matrix holds them before going on."
        )

    # --- column check ---------------------------------------------------
    missing = {k: v for k, v in cols.items() if k != "gene_symbol" and v not in obs.columns}
    if cols.get("gene_symbol") not in var.columns:
        missing["gene_symbol"] = cols.get("gene_symbol")
    if missing:
        out.append("\n## Missing columns\n")
        for k, v in missing.items():
            pool = var.columns if k == "gene_symbol" else obs.columns
            close = difflib.get_close_matches(str(v), list(pool), n=3, cutoff=0.4)
            out.append(f"- `{k}` -> `{v}` not found. Closest: {close or 'none'}")
        out.append("\nFix the names under `columns:` in config/config.yaml.")
    out.append("\n## All obs columns\n")
    out.append(", ".join(f"`{c}`" for c in obs.columns))

    c = {k: v for k, v in cols.items() if v in obs.columns}
    if missing and any(k in missing for k in ("donor", "disease")):
        report = "\n".join(out)
        print(report)
        sys.exit("Donor or disease column missing; fix config before continuing.")

    # --- samples and donors ---------------------------------------------
    key = [k for k in cfg["sample_key"] if k in obs.columns]
    obs["sample_id"] = make_sample_ids(obs, key).values
    samples = obs.groupby("sample_id", observed=True).agg(
        donor=(c["donor"], "first"),
        disease=(c["disease"], "first"),
        n_cells=(c["donor"], "size"),
        **({"n_libraries": (c["library"], "nunique")} if "library" in c else {}),
        **({"n_batches": (c["batch"], "nunique")} if "batch" in c else {}),
    )
    out.append("\n## Design\n")
    out.append(f"- Donors: **{obs[c['donor']].nunique()}**")
    out.append(f"- Samples (`{' x '.join(key)}`): **{len(samples)}**")
    per_donor = samples.groupby("donor", observed=True).size().value_counts().sort_index()
    out.append("\nSamples per donor:\n")
    out.append(md_table(per_donor.rename_axis("samples_per_donor").rename("donors")))

    donors = obs.groupby(c["donor"], observed=True).first()
    out.append("\n### Donors by disease\n")
    out.append(md_table(donors[c["disease"]].value_counts()))
    if "disease_state" in c:
        out.append("\n### Disease x disease_state (samples)\n")
        st = obs.groupby("sample_id", observed=True)[c["disease_state"]].first()
        out.append(md_table(pd.crosstab(samples["disease"], st)))

    # --- confounding ----------------------------------------------------
    out.append("\n## Confounding check (donors)\n")
    out.append(
        "Uneven rows here are covariates a classifier could exploit instead of "
        "disease. They must go into the DE model and be tested in validation."
    )
    for k in ("batch", "ancestry", "sex"):
        if k in c:
            # donor-level value; donors seen in several batches count as "mixed"
            per = obs.groupby(c["donor"], observed=True)[c[k]].agg(
                lambda s: s.iloc[0] if s.nunique() == 1 else "mixed"
            )
            tab = pd.crosstab(per, donors[c["disease"]], margins=True)
            out.append(f"\n### {k} (`{c[k]}`) x disease\n")
            out.append(md_table(tab))

    if "n_batches" in samples:
        multi = samples[samples["n_batches"] > 1]
        out.append(
            f"\nSamples spread over more than one batch (technical replicates): "
            f"**{len(multi)}**. Their libraries are summed into one pseudobulk "
            "profile per sample, and CV folds are grouped by donor so replicates "
            "never sit on both sides of a split."
        )

    # --- cells -----------------------------------------------------------
    out.append("\n## Cells per sample\n")
    out.append(md_table(samples["n_cells"].describe().round(0)))
    if "cell_type" in c:
        out.append("\n## Cell types (% of cells within each disease group)\n")
        ct = pd.crosstab(obs[c["cell_type"]], obs[c["disease"]], normalize="columns")
        out.append(md_table((ct * 100).round(2)))
        out.append(
            "\nCell-level percentages are dominated by samples with many cells; "
            "the proportion analysis in phase 2 is done per sample instead."
        )

    report = "\n".join(out) + "\n"
    rdir = Path(cfg["reports_dir"])
    rdir.mkdir(parents=True, exist_ok=True)
    (rdir / "01_data_audit.md").write_text(report)
    samples.to_csv(rdir / "01_samples.csv")
    print(report)
    print(f"Written: {rdir / '01_data_audit.md'} and {rdir / '01_samples.csv'}")


if __name__ == "__main__":
    main()
