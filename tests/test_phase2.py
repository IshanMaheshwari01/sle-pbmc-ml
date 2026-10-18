"""Phase 2 end to end on the mock: export -> R proportions -> R differential expression.

The mock plants two known effects, so the tests check the statistics recover
them and do not invent others:
* lupus donors have more classical monocytes and fewer CD4 T cells
* ten interferon genes are 4x higher in lupus monocytes, and only there

Skipped automatically when Rscript or the R packages are missing.
Run:  pytest -q tests/test_phase2.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mock_data import build_mock  # noqa: E402

if not shutil.which("Rscript"):
    pytest.skip("Rscript not installed", allow_module_level=True)
_check = subprocess.run(
    ["Rscript", "-e", "for (p in c('limma','edgeR','ggplot2','yaml','Matrix')) "
     "if (!requireNamespace(p, quietly=TRUE)) quit(status=1)"], capture_output=True)
if _check.returncode != 0:
    pytest.skip("R packages for phase 2 not installed", allow_module_level=True)

IFN = ["IFI27", "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1", "MX1", "OAS1", "IFI6", "IFIT3"]


def _run(cmd):
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    assert res.returncode == 0, f"{cmd[1]} failed:\n{res.stdout[-3000:]}\n{res.stderr[-5000:]}"
    return res.stdout


@pytest.fixture(scope="module")
def phase2(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("p2")
    (tmp / "config").mkdir()
    raw = tmp / "data/raw/mock.h5ad"
    raw.parent.mkdir(parents=True)
    build_mock(n_sle=26, n_healthy=22, seed=5).write_h5ad(raw)

    cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    cfg["data"]["raw_h5ad"] = str(raw)
    cfg["subsample"].update(max_cells_per_sample=40, out_h5ad=str(tmp / "data/processed/cells.h5ad"))
    cfg["pseudobulk"].update(out_h5ad=str(tmp / "data/processed/pb.h5ad"), split_by=["Processing_Cohort"])
    cfg["reports_dir"] = str(tmp / "results/reports")
    (tmp / "config/config.yaml").write_text(yaml.safe_dump(cfg))

    gmt = tmp / "mock_sets.gmt"
    lines = ["MOCK_INTERFERON\tplanted\t" + "\t".join(IFN)]
    for k in range(7):  # the cell-type identity blocks of the mock
        lines.append(f"MOCK_BLOCK_{k}\tidentity\t" + "\t".join(f"GENE{i}" for i in range(40 + 20 * k, 60 + 20 * k)))
    gmt.write_text("\n".join(lines) + "\n")
    p2 = yaml.safe_load((ROOT / "config/phase2.yaml").read_text())
    p2["cohort"]["min_ancestry_donors"] = 3
    p2["proportions"]["min_cells_per_profile"] = 20
    p2["de"].update(min_cells=5, min_profiles_per_group=5)
    p2["enrichment"].update(gmt=str(gmt), min_set_size=5)
    (tmp / "config/phase2.yaml").write_text(yaml.safe_dump(p2))

    _run([sys.executable, ROOT / "scripts/02_subsample.py", "--config", tmp / "config/config.yaml"])
    _run([sys.executable, ROOT / "scripts/06_export_pseudobulk.py", "--config", tmp / "config/config.yaml",
          "--phase2-config", tmp / "config/phase2.yaml"])
    data_dir = tmp / "data/processed/phase2"
    _run(["Rscript", ROOT / "R/07_proportions.R", data_dir])
    _run(["Rscript", ROOT / "R/08_differential_expression.R", data_dir])
    return tmp


def test_cohort_rules(phase2):
    prof = pd.read_csv(phase2 / "data/processed/phase2/profiles.csv")
    co = prof[prof["in_cohort"]]
    assert set(co["sex"]) == {"female"}
    assert co.groupby("donor")["sample_id"].nunique().max() == 1  # one sample per donor
    assert (co["batch"] != "mixed").all()
    assert (phase2 / "results/reports/06_cohort.md").exists()


def test_proportions_recover_planted_shift(phase2):
    res = pd.read_csv(phase2 / "results/reports/07_proportions.csv").set_index("cell_type")
    assert res.loc["classical monocyte", "odds_ratio"] > 1 and res.loc["classical monocyte", "fdr"] < 0.05
    assert res.loc["CD4-positive, alpha-beta T cell", "odds_ratio"] < 1
    assert res.loc["CD4-positive, alpha-beta T cell", "fdr"] < 0.05
    # no shift was planted for B cells
    assert res.loc["B cell", "fdr"] > 0.05
    assert (phase2 / "results/figures/07_proportions_forest.png").stat().st_size > 10_000


def test_de_finds_interferon_in_monocytes_only(phase2):
    de = phase2 / "results/reports/08_de"
    mono = pd.read_csv(de / "classical_monocyte_limma.csv").set_index("symbol")
    hits = mono.loc[mono.index.intersection(IFN)]
    assert len(hits) >= 8
    assert (hits["adj.P.Val"] < 0.05).mean() >= 0.8 and (hits["logFC"] > 1).mean() >= 0.8
    cd4 = pd.read_csv(de / "cd4_positive_alpha_beta_t_cell_limma.csv").set_index("symbol")
    cd4_ifn = cd4.loc[cd4.index.intersection(IFN)]
    assert (cd4_ifn["adj.P.Val"] < 0.05).mean() <= 0.2, "interferon genes called in CD4 T cells, where none was planted"


def test_deseq2_agrees(phase2):
    summ = pd.read_csv(phase2 / "results/reports/08_de_summary.csv").set_index("cell_type")
    assert summ.loc["classical monocyte", "deseq2_rho"] > 0.5
    assert summ.loc["classical monocyte", "pct_limma_de_also_deseq2"] >= 50


def test_pathway_test_finds_planted_set(phase2):
    cam = pd.read_csv(phase2 / "results/reports/08_pathways_camera.csv")
    row = cam[(cam["cell_type"] == "classical monocyte") & (cam["pathway"] == "MOCK_INTERFERON")]
    assert len(row) == 1 and row["Direction"].iloc[0] == "Up" and row["FDR"].iloc[0] < 0.05
    for f in ("08_de_counts.png", "08_pathways_heatmap.png", "08_limma_vs_deseq2.png"):
        assert (phase2 / "results/figures" / f).stat().st_size > 10_000
    report = (phase2 / "results/reports/08_differential_expression.md").read_text()
    assert "Interferon genes" in report
