"""End-to-end tests on a mock dataset with the real file's layout.

Run:  pytest -q
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mock_data import build_mock  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CAP, DEV_CAP = 60, 15


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("proj")
    src = build_mock(seed=1)
    raw_path = tmp / "data/raw/mock.h5ad"
    raw_path.parent.mkdir(parents=True)
    src.write_h5ad(raw_path)

    cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    cfg["data"]["raw_h5ad"] = str(raw_path)
    cfg["subsample"].update(max_cells_per_sample=CAP, block_size=37,  # odd size: tests block edges
                            out_h5ad=str(tmp / "sub.h5ad"))
    cfg["dev"].update(max_cells_per_sample=DEV_CAP, out_h5ad=str(tmp / "dev.h5ad"))
    cfg["pseudobulk"]["out_h5ad"] = str(tmp / "pb.h5ad")
    cfg["reports_dir"] = str(tmp / "reports")
    cfg_path = tmp / "config" / "config.yaml"
    cfg_path.parent.mkdir()
    cfg_path.write_text(yaml.safe_dump(cfg))

    for script, extra in [("01_inspect.py", []), ("02_subsample.py", []), ("02_subsample.py", ["--dev"])]:
        res = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / script), "--config", str(cfg_path), *extra],
            capture_output=True, text=True,
        )
        assert res.returncode == 0, f"{script} failed:\n{res.stdout}\n{res.stderr}"

    # second pseudobulk run with technical replicates summed instead of split
    cfg["pseudobulk"].update(split_by=[], out_h5ad=str(tmp / "pb_summed.h5ad"))
    cfg["subsample"]["out_h5ad"] = str(tmp / "sub2.h5ad")
    cfg_path2 = tmp / "config" / "summed.yaml"
    cfg_path2.write_text(yaml.safe_dump(cfg))
    res = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "02_subsample.py"), "--config", str(cfg_path2)],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr

    sub = ad.read_h5ad(tmp / "sub.h5ad")
    dev = ad.read_h5ad(tmp / "dev.h5ad")
    pb = ad.read_h5ad(tmp / "pb.h5ad")
    pb_summed = ad.read_h5ad(tmp / "pb_summed.h5ad")
    raw = ad.AnnData(X=src.raw.X, obs=src.obs.copy(), var=src.raw.var.copy())
    raw.obs["sample_id"] = raw.obs["donor_id"].astype(str) + "|" + raw.obs["disease_state"].astype(str)
    raw.obs["unit_id"] = raw.obs["sample_id"] + "|" + raw.obs["Processing_Cohort"].astype(str)
    return dict(tmp=tmp, raw=raw, sub=sub, dev=dev, pb=pb, pb_summed=pb_summed)


def test_every_sample_kept_and_capped(run):
    raw, sub, dev = run["raw"], run["sub"], run["dev"]
    all_samples = set(raw.obs["sample_id"])
    for obj, cap in [(sub, CAP), (dev, DEV_CAP)]:
        sizes = obj.obs["sample_id"].value_counts()
        assert set(sizes.index) == all_samples
        assert sizes.max() <= cap
        # a sample smaller than the cap keeps all its cells
        full = raw.obs["sample_id"].value_counts()
        small = full[full <= cap]
        assert (sizes.loc[small.index] == small).all()


def test_subsample_counts_identical_to_source(run):
    raw, sub = run["raw"], run["sub"]
    src = raw[sub.obs_names]
    assert np.array_equal(src.X.toarray(), sub.X.toarray())
    assert np.array_equal(sub.layers["counts"].toarray(), sub.X.toarray())
    assert list(sub.var_names) == list(raw.var_names)
    assert "feature_name" in sub.var


@pytest.mark.parametrize("which,unit", [("pb", "unit_id"), ("pb_summed", "sample_id")])
def test_pseudobulk_matches_bruteforce(run, which, unit):
    raw, pb = run["raw"], run[which]
    assert pb.X.sum() == pytest.approx(raw.X.sum())
    key = raw.obs[unit] + "||" + raw.obs["cell_type"].astype(str)
    dense = pd.DataFrame(raw.X.toarray(), index=raw.obs_names)
    expected = dense.groupby(key.values).sum()
    got = pd.DataFrame(pb.X.toarray(), index=pb.obs_names)
    assert set(expected.index) == set(got.index)
    assert np.allclose(expected.loc[got.index].to_numpy(), got.to_numpy())
    assert (pb.obs["n_cells"] == key.value_counts().loc[pb.obs_names]).all()


def test_pseudobulk_split_by_batch(run):
    pb = run["pb"]
    # every profile now comes from a single batch
    assert (pb.obs["batch"] != "mixed").all()
    assert (pb.obs["n_batches"] == 1).all()
    # cross-batch replicates of one sample appear as separate profiles
    per = pb.obs.groupby(["sample_id", "cell_type"], observed=True)["profile_id"].nunique()
    assert (per > 1).any()
    # one donor, two time points -> two samples of the same donor
    multi = pb.obs.groupby("donor", observed=True)["sample_id"].nunique()
    assert (multi > 1).any()
    assert pb.obs["low_cell_count"].dtype == bool


def test_pseudobulk_summed_flags_mixed_batches(run):
    pb = run["pb_summed"]
    assert (pb.obs["batch"] == "mixed").any()
    assert (pb.obs.loc[pb.obs["batch"] == "mixed", "n_batches"] > 1).all()


def test_audit_report_written(run):
    report = (run["tmp"] / "reports" / "01_data_audit.md").read_text()
    assert "Confounding check" in report
    assert "all whole numbers" in report
    assert "Missing columns" not in report
