"""Load config/config.yaml and resolve its paths."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "config.yaml"

# Any key ending in one of these is a path. Relative paths are resolved
# against the project root (the folder that contains config/), so scripts work
# from any working directory.
_PATH_SUFFIXES = ("_h5ad", "_dir", "_csv")


def _resolve(node: Any, root: Path) -> Any:
    if isinstance(node, dict):
        out = {}
        for k, v in node.items():
            if isinstance(v, str) and k.endswith(_PATH_SUFFIXES):
                p = Path(v)
                out[k] = p if p.is_absolute() else root / p
            else:
                out[k] = _resolve(v, root)
        return out
    return node


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    root = path.resolve().parent.parent
    cfg = _resolve(cfg, root)
    cfg["_root"] = root
    return cfg
