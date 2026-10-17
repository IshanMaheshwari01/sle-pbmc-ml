"""Read a very large .h5ad without loading its count matrix into memory.

The CELLxGENE file for this dataset is about 12 GB. AnnData's backed mode does
not handle `.raw` (where CELLxGENE keeps the raw counts) well, so this module
reads the HDF5 groups directly:

* `obs` and `var` are small and are read whole.
* The count matrix is CSR (one row per cell), so it is read in contiguous
  blocks of rows. Each block costs one sequential read, which is far faster
  than thousands of single-row reads from a compressed file.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:  # anndata >= 0.11
    from anndata.io import read_elem
except ImportError:  # pragma: no cover - older anndata
    from anndata.experimental import read_elem


@dataclass
class CountsLocation:
    """Where the raw counts live inside the file."""

    matrix_path: str  # e.g. "raw/X"
    var_path: str  # e.g. "raw/var"
    shape: tuple[int, int]


def read_obs(path) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


def locate_counts(path) -> CountsLocation:
    """Find the raw-count matrix.

    CELLxGENE stores normalised values in X and raw counts in raw/X. Some
    re-processed copies keep counts in layers/counts or directly in X, so
    those are checked too, in that order.
    """
    with h5py.File(path, "r") as f:
        candidates = [
            ("raw/X", "raw/var"),
            ("layers/counts", "var"),
            ("X", "var"),
        ]
        for mat, var in candidates:
            if mat in f and var in f:
                node = f[mat]
                if not isinstance(node, h5py.Group):
                    raise ValueError(
                        f"{mat} is a dense array; expected a sparse CSR matrix."
                    )
                enc = node.attrs.get("encoding-type", "")
                enc = enc.decode() if isinstance(enc, bytes) else enc
                if enc != "csr_matrix":
                    raise ValueError(
                        f"{mat} is stored as {enc!r}; streaming needs csr_matrix "
                        "(one row per cell)."
                    )
                shape = tuple(int(x) for x in node.attrs["shape"])
                return CountsLocation(mat, var, shape)
    raise ValueError(f"No sparse count matrix found in {path}.")


def read_var(path, loc: CountsLocation) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        var = read_elem(f[loc.var_path])
        # raw/var sometimes lacks the annotation columns that the main var has
        # (CELLxGENE keeps feature_name in both, but be defensive).
        if loc.var_path != "var" and "var" in f:
            main_var = read_elem(f["var"])
            extra = [c for c in main_var.columns if c not in var.columns]
            if extra:
                var = var.join(main_var[extra], how="left")
    return var


def iter_row_blocks(
    path, loc: CountsLocation, block_size: int
) -> Iterator[tuple[int, sp.csr_matrix]]:
    """Yield (first_row, csr_block) covering every row in order."""
    n_rows, n_cols = loc.shape
    with h5py.File(path, "r") as f:
        g = f[loc.matrix_path]
        indptr = g["indptr"][:]
        data_ds, idx_ds = g["data"], g["indices"]
        for start in range(0, n_rows, block_size):
            stop = min(start + block_size, n_rows)
            lo, hi = int(indptr[start]), int(indptr[stop])
            block = sp.csr_matrix(
                (data_ds[lo:hi], idx_ds[lo:hi], indptr[start : stop + 1] - lo),
                shape=(stop - start, n_cols),
            )
            yield start, block


def check_integer_counts(path, loc: CountsLocation, n_values: int = 200_000) -> dict:
    """Look at the first stored values and report whether they are whole numbers."""
    with h5py.File(path, "r") as f:
        ds = f[loc.matrix_path]["data"]
        disk_dtype = ds.dtype
        vals = np.asarray(ds[:n_values], dtype=np.float64)
    is_int = bool(np.all(np.isclose(vals, np.round(vals))))
    return {
        "checked_values": int(vals.size),
        "all_integer": is_int,
        "min": float(vals.min()) if vals.size else None,
        "max": float(vals.max()) if vals.size else None,
        "dtype_on_disk": str(disk_dtype),
    }
