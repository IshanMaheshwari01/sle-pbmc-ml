"""Phase 3 on the mock: the classifier, its guards against leakage, and its stress tests.

The mock plants an interferon signal in lupus monocytes and a separate batch
effect, so the tests check that:
* the main model finds lupus within every batch, not only when pooled
* gene selection happens inside folds and no donor is on both sides of a split
* a shuffled-label run scores near chance
* the batch check detects the planted batch effect among healthy donors
Run:  pytest -q tests/test_phase3.py
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
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mock_data import build_mock  # noqa: E402
from sle_pbmc.phase3 import (  # noqa: E402
    BlockTopVariance,
    build_dataset,
    model_specs,
    read_phase2,
    run_cv,
)

pytest.importorskip("sklearn")


def _run(cmd):
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    assert res.returncode == 0, f"{cmd[1]} failed:\n{res.stdout[-3000:]}\n{res.stderr[-5000:]}"
    return res.stdout


@pytest.fixture(scope="module")
def phase3(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("p3")
    (tmp / "config").mkdir()
    raw = tmp / "data/raw/mock.h5ad"
    raw.parent.mkdir(parents=True)
    build_mock(n_sle=40, n_healthy=34, seed=11).write_h5ad(raw)

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
    p3["features"]["min_cells"] = 3
    p3["features"]["top_genes_per_block"] = 100
    p3["cv"]["repeats"] = 2
    p3["evaluation"].update(permutations=10, bootstrap=200)
    (tmp / "config/phase3.yaml").write_text(yaml.safe_dump(p3))

    _run([sys.executable, ROOT / "scripts/02_subsample.py", "--config", tmp / "config/config.yaml"])
    _run([sys.executable, ROOT / "scripts/06_export_pseudobulk.py", "--config", tmp / "config/config.yaml",
          "--phase2-config", tmp / "config/phase2.yaml"])
    _run([sys.executable, ROOT / "scripts/09_classifier.py", "--config", tmp / "config/config.yaml",
          "--phase3-config", tmp / "config/phase3.yaml"])
    counts, profiles, genes = read_phase2(tmp / "data/processed/phase2")
    ds = build_dataset(counts, profiles, genes, p3)
    return {"tmp": tmp, "ds": ds, "p3": p3}


def test_dataset_shape(phase3):
    ds = phase3["ds"]
    assert len(ds.y) == len(ds.meta) == len(ds.X)
    assert ds.meta.groupby("donor")["sample_id"].nunique().max() == 1  # one sample per donor
    assert "ifn|classical monocyte" in ds.X.columns
    assert all(c.startswith("cov|") for c in ds.blocks["covariates"])
    full = next(s for s in model_specs(ds) if s.name.startswith("L1"))
    assert not any(c.startswith("cov|") or c.startswith("ifn|") for c in full.columns)


def test_gene_selection_uses_training_rows_only():
    X = pd.DataFrame({"expr|a|g1": [0, 0, 0, 10.0], "expr|a|g2": [1, 2, 3, 4.0], "prop|x": [1, 2, 3, 4.0]})
    sel = BlockTopVariance({"a expression": ["expr|a|g1", "expr|a|g2"]}, n=1)
    assert list(sel.fit(X.iloc[:3]).columns_) == ["expr|a|g2", "prop|x"]  # g1 is constant in training rows
    assert list(sel.fit(X).columns_) == ["expr|a|g1", "prop|x"]


def test_no_donor_on_both_sides(phase3):
    from sklearn.model_selection import StratifiedGroupKFold

    ds = phase3["ds"]
    g = ds.meta["donor"].to_numpy()
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=0).split(ds.X, ds.y, g):
        assert not set(g[tr]) & set(g[te])


def test_report_and_within_batch_accuracy(phase3):
    perf = pd.read_csv(phase3["tmp"] / "results/reports/09_performance.csv")
    main = perf[perf["model"] == "L1 logistic (all features)"].set_index("setting")
    assert main.loc["All donors (pooled)", "AUROC"] > 0.9
    batch_rows = [s for s in main.index if s.startswith("Batch")]
    assert batch_rows and (main.loc[batch_rows, "AUROC"] > 0.85).all()
    cov = perf[perf["model"].str.startswith("Covariates")]
    assert not cov["setting"].str.startswith("Batch").any()  # not scored within a batch
    report = (phase3["tmp"] / "results/reports/09_classifier.md").read_text()
    for s in ("Permutation test", "How much batch is written", "without a high interferon"):
        assert s in report
    for f in ("09_auroc_by_setting.png", "09_roc_within_batch.png", "09_feature_block_importance.png",
              "09_predictions_by_ifn_status.png"):
        assert (phase3["tmp"] / "results/figures" / f).stat().st_size > 10_000


def test_shuffled_labels_score_near_chance(phase3):
    ds, p3 = phase3["ds"], phase3["p3"]
    spec = next(s for s in model_specs(ds) if s.name.startswith("L1"))
    rng = np.random.default_rng(1)
    donors = ds.meta["donor"].to_numpy()
    dy = pd.Series(ds.y, index=donors).groupby(level=0).first()
    aucs = []
    from sklearn.metrics import roc_auc_score

    for i in range(5):
        perm = pd.Series(rng.permutation(dy.to_numpy()), index=dy.index).loc[donors].to_numpy()
        r = run_cv(spec, ds, p3, seed=50 + i, repeats=1, y=perm, C_fixed=0.1, verbose=False)
        aucs.append(roc_auc_score(perm, r.preds[0]))
    assert 0.3 < np.mean(aucs) < 0.65


def test_batch_check_detects_planted_batch_effect(phase3):
    report = (phase3["tmp"] / "results/reports/09_classifier.md").read_text()
    line = next(l for l in report.splitlines() if "separate batch" in l)
    auc = float(line.split("AUROC **")[1].split("**")[0])
    assert auc > 0.85
