"""Phase 1 end to end on the mock dataset: QC -> integration -> CellTypist check.

Skipped automatically when the phase-1 packages are not installed, so the
phase-0 environment still passes `make test`.
Run:  pytest -q tests/test_phase1.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

for mod in ("scanpy", "scvi", "harmonypy", "celltypist", "skmisc"):
    pytest.importorskip(mod)

import anndata as ad  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mock_data import build_mock  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

# mock author labels -> CellTypist-style names, to train a tiny local model
CT_NAMES = {
    "CD4-positive, alpha-beta T cell": "Tcm/Naive helper T cells",
    "CD8-positive, alpha-beta T cell": "Tem/Temra cytotoxic T cells",
    "classical monocyte": "Classical monocytes",
    "non-classical monocyte": "Non-classical monocytes",
    "B cell": "Naive B cells",
    "natural killer cell": "CD16+ NK cells",
    "plasmacytoid dendritic cell": "pDC",
}


def _run(script, *args):
    res = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
                         capture_output=True, text=True)
    assert res.returncode == 0, f"{script} failed:\n{res.stdout[-3000:]}\n{res.stderr[-5000:]}"
    return res.stdout


@pytest.fixture(scope="module")
def phase1(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("p1")
    (tmp / "config").mkdir()
    raw = tmp / "data/raw/mock.h5ad"
    raw.parent.mkdir(parents=True)
    build_mock(n_sle=16, n_healthy=12, seed=3).write_h5ad(raw)

    cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    cfg["data"]["raw_h5ad"] = str(raw)
    cfg["subsample"].update(max_cells_per_sample=120, out_h5ad=str(tmp / "data/processed/cells.h5ad"))
    cfg["dev"]["out_h5ad"] = str(tmp / "data/processed/dev.h5ad")
    cfg["pseudobulk"]["out_h5ad"] = str(tmp / "data/processed/pb.h5ad")
    cfg["reports_dir"] = str(tmp / "results/reports")
    cfg_path = tmp / "config/config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    p1 = yaml.safe_load((ROOT / "config/phase1.yaml").read_text())
    p1["qc"]["min_genes"] = 20  # the mock has only 300 genes
    p1["hvg"]["n_top_genes"] = 200
    p1["pca"]["n_comps"] = 15
    p1["scvi"].update(n_latent=5, n_hidden=32, n_layers=1, batch_size=128)
    p1["benchmark"]["n_cells"] = 1500
    p1["disease_signal"].update(batch="auto", min_cells_per_group=10)
    p1["celltypist"]["chunk_size"] = 1000  # forces several chunks
    p1_path = tmp / "config/phase1.yaml"
    p1_path.write_text(yaml.safe_dump(p1))

    common = ["--config", cfg_path, "--phase1-config", p1_path]
    _run("02_subsample.py", "--config", cfg_path)
    _run("03_qc.py", *common)
    _run("04_integrate.py", *common, "--max-epochs", 3)

    # train a tiny CellTypist model on the mock so the annotation step runs offline
    import celltypist

    integ = ad.read_h5ad(tmp / "data/processed/cells_integrated.h5ad")
    q = ad.AnnData(X=integ.X, obs=pd.DataFrame(
        {"label": integ.obs["cell_type"].astype(str).map(CT_NAMES).values}, index=integ.obs_names))
    q.var_names = integ.var["feature_name"].astype(str).values
    model = celltypist.train(q, labels="label", n_jobs=1, max_iter=50, check_expression=False)
    model_path = tmp / "mock_model.pkl"
    model.write(str(model_path))
    _run("05_annotate.py", *common, "--model-path", model_path)
    return tmp


def test_qc_outputs(phase1):
    qc = ad.read_h5ad(phase1 / "data/processed/cells_qc.h5ad")
    assert qc.obs["qc_pass"].all()
    assert {"pct_counts_mt", "n_genes_by_counts"} <= set(qc.obs.columns)
    report = (phase1 / "results/reports/03_qc.md").read_text()
    assert "Flagged" in report and "Mitochondrial genes found: 5" in report
    assert (phase1 / "results/figures/03_qc_by_batch.png").stat().st_size > 10_000


def test_integration_outputs(phase1):
    a = ad.read_h5ad(phase1 / "data/processed/cells_integrated.h5ad")
    for key in ("X_pca", "X_harmony", "X_scVI", "X_umap_pca", "X_umap_harmony", "X_umap_scVI"):
        assert key in a.obsm and a.obsm[key].shape[0] == a.n_obs
        assert np.isfinite(a.obsm[key]).all()
    assert "leiden_scvi" in a.obs and "ifn_score" in a.obs
    assert "counts" in a.layers
    assert (phase1 / "results/models/cells/scvi/model.pt").exists()
    for f in ("04_umap_integration.png", "04_umap_cell_types.png", "04_ifn_score.png"):
        assert (phase1 / "results/figures" / f).stat().st_size > 10_000


def test_ifn_signal_detected(phase1):
    """The mock puts a 4x interferon signal in lupus monocytes; the score must find it."""
    report = (phase1 / "results/reports/04_integration.md").read_text()
    assert "Interferon score" in report
    a = ad.read_h5ad(phase1 / "data/processed/cells_integrated.h5ad")
    mono = a.obs["cell_type"].astype(str) == "classical monocyte"
    sle = a.obs["disease"].astype(str) == "systemic lupus erythematosus"
    assert a.obs.loc[mono & sle, "ifn_score"].median() > a.obs.loc[mono & ~sle, "ifn_score"].median()


def test_disease_signal_table(phase1):
    report = (phase1 / "results/reports/04_integration.md").read_text()
    assert "Is disease signal still there?" in report
    assert "Unintegrated (PCA)" in report and "scVI" in report


def test_annotation_outputs(phase1):
    pred = pd.read_csv(phase1 / "results/reports/05_celltypist_predictions.csv", index_col=0)
    assert len(pred) > 0
    agree = (pred["celltypist_lineage"] == pred["author_lineage"]).mean()
    assert agree > 0.8, f"lineage agreement only {agree:.2f} on easy mock data"
    assert (phase1 / "results/figures/05_celltypist_vs_author.png").stat().st_size > 10_000
