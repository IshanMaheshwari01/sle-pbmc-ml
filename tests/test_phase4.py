"""Phase 4 pool check on the mock.

* pooled mock (donors mixed across libraries): within-library AUROC computed,
  libraries reported as mixed, README section filled in
* the within-library AUROC and library-mate metrics on hand-made examples
Run:  pytest -q tests/test_phase4.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "scripts"))
from mock_data import build_mock  # noqa: E402

pc = pytest.importorskip("importlib").import_module("10_pool_check")


def _run(cmd):
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    assert res.returncode == 0, f"{cmd[1]} failed:\n{res.stdout[-3000:]}\n{res.stderr[-5000:]}"
    return res.stdout


def test_within_library_auc_by_hand():
    units = pd.DataFrame({"disease_group": ["SLE", "SLE", "healthy", "healthy"]}, index=list("abcd"))
    libs = pd.Series({"a": ["L1"], "b": ["L2"], "c": ["L1"], "d": ["L2"]})
    # within L1: a (0.9) vs c (0.2) -> concordant; within L2: b (0.1) vs d (0.8) -> discordant
    p = pd.Series({"a": 0.9, "b": 0.1, "c": 0.2, "d": 0.8})
    auc, n, per_lib = pc.within_library_auc(units, p, libs)
    assert n == 2 and auc == pytest.approx(0.5) and len(per_lib) == 2
    # libraries holding a single group contribute no pairs
    libs2 = pd.Series({"a": ["L1"], "b": ["L1"], "c": ["L2"], "d": ["L2"]})
    assert np.isnan(pc.within_library_auc(units, p, libs2)[0])
    # separated libraries: nobody has a partner, Cramer's V = 1; mixed: everyone does, V = 0
    assert pc.library_mixing(units, libs2) == (pytest.approx(0.0), pytest.approx(1.0))
    assert pc.library_mixing(units, libs) == (pytest.approx(100.0), pytest.approx(0.0))
    assert "largely disappears" in pc.verdict(0.6) and "holds" in pc.verdict(0.9)


@pytest.fixture(scope="module")
def phase4(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("p4")
    (tmp / "config").mkdir()
    raw = tmp / "data/raw/mock.h5ad"
    raw.parent.mkdir(parents=True)
    build_mock(n_sle=30, n_healthy=26, seed=21, pooled=True).write_h5ad(raw)
    cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    cfg["data"]["raw_h5ad"] = str(raw)
    cfg["subsample"].update(max_cells_per_sample=20, out_h5ad=str(tmp / "data/processed/cells.h5ad"))
    cfg["pseudobulk"].update(out_h5ad=str(tmp / "data/processed/pb.h5ad"), split_by=["Processing_Cohort"])
    cfg["reports_dir"] = str(tmp / "results/reports")
    (tmp / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    p2 = yaml.safe_load((ROOT / "config/phase2.yaml").read_text())
    p2["cohort"]["min_ancestry_donors"] = 3
    (tmp / "config/phase2.yaml").write_text(yaml.safe_dump(p2))
    p3 = yaml.safe_load((ROOT / "config/phase3.yaml").read_text())
    p3["features"].update(min_cells=3, top_genes_per_block=100)
    (tmp / "config/phase3.yaml").write_text(yaml.safe_dump(p3))
    readme = tmp / "README.md"
    readme.write_text("# x\n\n<!-- POOL_CHECK_START -->\nPLACEHOLDER_TEXT\n<!-- POOL_CHECK_END -->\n\nafter\n")
    (tmp / "SUMMARY.md").write_text("- <!-- POOL_SENTENCE_START -->PLACEHOLDER_TEXT<!-- POOL_SENTENCE_END -->\n")

    _run([sys.executable, ROOT / "scripts/02_subsample.py", "--config", tmp / "config/config.yaml"])
    _run([sys.executable, ROOT / "scripts/06_export_pseudobulk.py", "--config", tmp / "config/config.yaml",
          "--phase2-config", tmp / "config/phase2.yaml"])
    _run([sys.executable, ROOT / "scripts/10_pool_check.py", "--config", tmp / "config/config.yaml",
          "--phase3-config", tmp / "config/phase3.yaml", "--repeats", "1", "--readme", readme,
          "--summary", tmp / "SUMMARY.md"])
    return tmp


def test_pool_check_outputs(phase4):
    tab = pd.read_csv(phase4 / "results/reports/10_pool_check.csv").set_index("setting")
    allb = tab.loc["All batches"]
    assert allb["libraries with both groups"] >= 4
    assert allb["within-library AUROC (model)"] > 0.85         # planted IFN signal survives within libraries
    assert allb["% units with an opposite-group partner"] == 100   # mock libraries are mixed by design
    assert allb["Cramer's V (library x diagnosis)"] < 0.3
    pred = pd.read_csv(phase4 / "results/reports/10_unit_predictions.csv")
    assert {"p_L1 logistic (all features)", "libraries"} <= set(pred.columns)


def test_readme_section_filled(phase4):
    txt = (phase4 / "README.md").read_text()
    assert "PLACEHOLDER_TEXT" not in txt and "Within-library AUROC" in txt and txt.rstrip().endswith("after")
    assert txt.count("<!-- POOL_CHECK_START -->") == 1
    s = (phase4 / "SUMMARY.md").read_text()
    assert "PLACEHOLDER_TEXT" not in s and "Same batch, same 10x library" in s and "holds" in s
