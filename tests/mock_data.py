"""Build a small h5ad with the same layout as the CELLxGENE lupus file.

Lets you run the whole pipeline (and the tests) in seconds, before or without
the 12 GB download:  python tests/mock_data.py data/raw/mock.h5ad

Layout mirrored from the real file:
* X holds log-normalised values; raw/X holds integer counts (CSR)
* var index = Ensembl IDs, var["feature_name"] = gene symbols
* obs has donor_id, disease, disease_state, Processing_Cohort, library_uuid,
  self_reported_ethnicity, sex, cell_type, author_cell_type, plus extras
* some lupus donors have two time points (flare -> treated)
* some samples were sequenced in two libraries/batches (technical replicates)
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

CELL_TYPES = [
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "classical monocyte",
    "non-classical monocyte",
    "B cell",
    "natural killer cell",
    "plasmacytoid dendritic cell",
]
CT_PROBS = np.array([0.38, 0.14, 0.2, 0.05, 0.12, 0.1, 0.01])


def build_mock(n_sle: int = 12, n_healthy: int = 8, n_genes: int = 300, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    rows = []
    donors = [(f"SLE_{i:03d}", "systemic lupus erythematosus") for i in range(n_sle)] + [
        (f"HC_{i:03d}", "normal") for i in range(n_healthy)
    ]
    lib_counter = 0
    for d_idx, (donor, disease) in enumerate(donors):
        ancestry = rng.choice(["European", "Asian"])
        sex = rng.choice(["female", "male"], p=[0.85, 0.15])
        if disease == "normal":
            states = ["na"]
        elif d_idx % 4 == 0:
            states = ["flare", "treated"]  # same donor, two time points
        else:
            states = ["managed"]
        for state in states:
            # every 3rd sample is a technical replicate across two batches
            n_libs = 2 if (d_idx % 3 == 0) else 1
            for lib in range(n_libs):
                lib_counter += 1
                batch = str(1 + (d_idx + lib) % 4)
                n = int(rng.integers(30, 160))
                ct = rng.choice(CELL_TYPES, size=n, p=CT_PROBS)
                for c in ct:
                    rows.append(
                        dict(
                            donor_id=donor,
                            ind_cov=donor,
                            disease=disease,
                            disease_state=state,
                            Processing_Cohort=batch,
                            library_uuid=f"lib{lib_counter:04d}",
                            self_reported_ethnicity=ancestry,
                            sex=sex,
                            cell_type=c,
                            author_cell_type=c.split(",")[0][:6],
                            is_primary_data=True,
                            assay="10x 3' v2",
                        )
                    )
    obs = pd.DataFrame(rows)
    # shuffle row order so samples are interleaved, as in a pooled 10x design
    obs = obs.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    obs.index = [f"cell{i:06d}" for i in range(len(obs))]
    for col in obs.columns:
        if obs[col].dtype == object:
            obs[col] = obs[col].astype("category")

    n = len(obs)
    base = rng.gamma(0.4, 1.0, size=n_genes)
    lam = np.outer(rng.gamma(5, 0.4, size=n), base)
    # interferon-like signal: first 10 genes up in lupus monocytes
    sle_mono = (obs["disease"] == "systemic lupus erythematosus").to_numpy() & obs[
        "cell_type"
    ].astype(str).str.contains("monocyte").to_numpy()
    lam[np.ix_(sle_mono, np.arange(10))] *= 4
    counts = sp.csr_matrix(rng.poisson(lam).astype(np.float32))

    var = pd.DataFrame(
        {"feature_name": [f"GENE{i}" for i in range(n_genes)], "feature_biotype": "gene"},
        index=[f"ENSG{i:011d}" for i in range(n_genes)],
    )
    raw = ad.AnnData(X=counts, obs=obs, var=var)
    norm = counts.multiply(1e4 / np.maximum(counts.sum(1), 1)).tocsr()
    norm.data = np.log1p(norm.data)
    adata = ad.AnnData(X=norm.astype(np.float32), obs=obs, var=var)
    adata.raw = raw
    adata.uns["title"] = "MOCK of multiplexed scRNA-seq of PBMCs from lupus samples"
    return adata


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "data/raw/mock.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    a = build_mock()
    a.write_h5ad(out)
    print(f"Wrote mock dataset {out}: {a.n_obs} cells x {a.n_vars} genes")
