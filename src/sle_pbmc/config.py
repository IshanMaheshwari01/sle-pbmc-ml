"""Load config/config.yaml and resolve its paths."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "config.yaml"

# Keys whose values are file or directory paths. Relative paths are resolved
# against the project root (the folder that contains config/), so scripts work
# from any working directory.
_PATH_KEYS = {
    ("data", "raw_h5ad"),
    ("subsample", "out_h5ad"),
    ("dev", "out_h5ad"),
    ("pseudobulk", "out_h5ad"),
    ("reports_dir",),
}


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    root = path.resolve().parent.parent
    for keys in _PATH_KEYS:
        node = cfg
        for k in keys[:-1]:
            node = node.get(k, {})
        if keys[-1] in node and node[keys[-1]] is not None:
            p = Path(node[keys[-1]])
            node[keys[-1]] = p if p.is_absolute() else root / p
    cfg["_root"] = root
    return cfg
