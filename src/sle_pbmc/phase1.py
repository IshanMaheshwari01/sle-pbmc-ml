"""Shared helpers for phase 1 (QC, integration, annotation)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PHASE1_CONFIG = Path(__file__).resolve().parents[2] / "config" / "phase1.yaml"


def load_phase1(path: str | Path | None = None) -> dict:
    with open(Path(path) if path else PHASE1_CONFIG) as fh:
        return yaml.safe_load(fh)


def phase1_paths(cfg: dict, dev: bool = False) -> dict[str, Path]:
    """All phase-1 paths, derived from config.yaml so mock and real runs never collide."""
    base = Path(cfg["dev"]["out_h5ad"] if dev else cfg["subsample"]["out_h5ad"])
    reports = Path(cfg["reports_dir"])
    parts = list(reports.parts)
    idx = len(parts) - 1 - parts[::-1].index("reports") if "reports" in parts else None
    figures = Path(*parts[:idx], "figures", *parts[idx + 1 :]) if idx is not None else reports / "figures"
    if dev:
        figures = figures / "dev"
        reports = reports / "dev"
    return {
        "input": base,
        "qc": base.with_name(base.stem + "_qc.h5ad"),
        "integrated": base.with_name(base.stem + "_integrated.h5ad"),
        "reports": reports,
        "figures": figures,
        "models": Path(cfg["_root"]) / "results" / "models" / base.stem,
        "celltypist_csv": reports / "05_celltypist_predictions.csv",
    }


def mad_outlier(x: pd.Series, n_mads: float, side: str = "both") -> pd.Series:
    """True where x is more than n_mads MADs from its median."""
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826  # consistent with SD for normal data
    if mad == 0:
        return pd.Series(False, index=x.index)
    z = (x - med) / mad
    if side == "upper":
        return z > n_mads
    if side == "lower":
        return z < -n_mads
    return z.abs() > n_mads


def stratified_subset(obs: pd.DataFrame, keys: list[str], n: int, seed: int) -> np.ndarray:
    """Row positions for ~n cells, sampled proportionally within each stratum."""
    if n >= len(obs):
        return np.arange(len(obs))
    rng = np.random.default_rng(seed)
    frac = n / len(obs)
    pos = np.arange(len(obs))
    out = []
    for _, rows in pd.Series(pos).groupby([obs[k].astype(str).values for k in keys]):
        take = max(1, int(round(len(rows) * frac)))
        out.append(rng.choice(rows.to_numpy(), size=min(take, len(rows)), replace=False))
    return np.sort(np.concatenate(out))


def gene_symbols(var: pd.DataFrame, col: str) -> pd.Index:
    return pd.Index(var[col].astype(str)) if col in var.columns else var.index


def write_report(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(lines) + "\n"
    path.write_text(text)
    print(text)
    print(f"Written: {path}")
