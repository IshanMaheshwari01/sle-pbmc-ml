"""Phase 3: features, models and cross-validation for the donor-level classifier.

Unit of analysis: one sample in one processing batch (a "unit"). Most donors
have one unit; donors sequenced in two batches have two. Every split is
grouped by donor, so no person is ever on both sides of a split.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
@dataclass
class Dataset:
    X: pd.DataFrame            # units x features (may contain NaN)
    meta: pd.DataFrame         # units: donor, disease_group, batch, ancestry, ...
    y: np.ndarray              # 1 = SLE
    blocks: dict[str, list[str]] = field(default_factory=dict)  # block name -> columns


def read_phase2(dir_: Path):
    dir_ = Path(dir_)
    with gzip.open(dir_ / "counts.mtx.gz", "rb") as fh:
        counts = sp.csc_matrix(scipy.io.mmread(fh))  # genes x profiles
    profiles = pd.read_csv(dir_ / "profiles.csv", dtype={"batch": str})
    genes = pd.read_csv(dir_ / "genes.csv")
    profiles["in_cohort"] = profiles["in_cohort"].astype(str).str.lower() == "true"
    return counts, profiles, genes


def build_dataset(counts, profiles: pd.DataFrame, genes: pd.DataFrame, p3: dict) -> Dataset:
    f = p3["features"]
    co = profiles[profiles["in_cohort"]].copy()
    meta = (co.groupby("profile_id", sort=True)
              .agg(donor=("donor", "first"), sample_id=("sample_id", "first"),
                   disease_group=("disease_group", "first"), batch=("batch", "first"),
                   ancestry=("ancestry", "first"), disease_state=("disease_state", "first"),
                   n_cells=("n_cells", "sum")))
    units = meta.index
    blocks: dict[str, list[str]] = {}
    parts = []

    # cell-type composition, centred log-ratio
    N = co.pivot_table(index="profile_id", columns="cell_type", values="n_cells", aggfunc="sum", fill_value=0)
    N = N.reindex(units).fillna(0)
    logp = np.log(N + f["pseudocount_clr"])
    clr = logp.sub(logp.mean(axis=1), axis=0)
    clr.columns = [f"prop|{c}" for c in clr.columns]
    parts.append(clr)
    blocks["cell-type proportions"] = list(clr.columns)

    # covariates (only for the covariate baseline)
    cov = pd.get_dummies(meta[["batch", "ancestry"]].astype(str), prefix=["cov|batch", "cov|ancestry"], dtype=float)
    parts.append(cov)
    blocks["covariates"] = list(cov.columns)

    symbols = genes["symbol"].astype(str).to_numpy()
    ifn_cols = []
    for ct in f["cell_types"]:
        sub = co[(co["cell_type"] == ct) & (co["n_cells"] >= f["min_cells"])]
        if sub.empty:
            continue
        idx = sub.index.to_numpy()
        C = counts[:, idx].T.tocsr().astype(np.float64)        # units x genes
        lib = np.asarray(C.sum(axis=1)).ravel()
        lib[lib == 0] = 1
        cpm = sp.diags(1e6 / lib) @ C
        detect = np.asarray((cpm >= 1).mean(axis=0)).ravel()
        keep = np.flatnonzero(detect >= f["min_detect_frac"])
        logcpm = np.log2(cpm[:, keep].toarray() + 1).astype(np.float32)
        names = symbols[keep]
        # unique column names even if a symbol repeats
        seen: dict[str, int] = {}
        cols = []
        for s in names:
            k = seen.get(s, 0)
            seen[s] = k + 1
            cols.append(f"expr|{ct}|{s}" + (f"#{k}" if k else ""))
        E = pd.DataFrame(logcpm, index=sub["profile_id"].to_numpy(), columns=cols).reindex(units)
        parts.append(E)
        blocks[f"{ct} expression"] = cols
        ifn_present = [i for i, s in enumerate(names) if s in set(f["ifn_genes"])]
        if ifn_present:
            ifn = pd.Series(logcpm[:, ifn_present].mean(axis=1), index=sub["profile_id"].to_numpy()).reindex(units)
            parts.append(ifn.rename(f"ifn|{ct}").to_frame())
            ifn_cols.append(f"ifn|{ct}")
    blocks["interferon scores"] = ifn_cols

    X = pd.concat(parts, axis=1)
    y = (meta["disease_group"] == "SLE").astype(int).to_numpy()
    return Dataset(X=X, meta=meta, y=y, blocks=blocks)


# --------------------------------------------------------------------------
# pipeline pieces
# --------------------------------------------------------------------------
class BlockTopVariance(BaseEstimator, TransformerMixin):
    """Keep the `n` most variable columns of each expression block, chosen on
    the training rows only; other columns pass through."""

    def __init__(self, blocks: dict | None = None, n: int = 1000):
        self.blocks = blocks
        self.n = n

    def fit(self, X, y=None):
        cols = np.asarray(X.columns)
        keep = np.ones(len(cols), dtype=bool)
        var = np.nanvar(X.to_numpy(dtype=np.float64), axis=0)
        var = np.nan_to_num(var, nan=-1.0)
        pos = {c: i for i, c in enumerate(cols)}
        for name, bcols in (self.blocks or {}).items():
            if not name.endswith("expression"):
                continue
            idx = np.array([pos[c] for c in bcols if c in pos])
            if len(idx) > self.n:
                keep[idx] = False
                keep[idx[np.argsort(var[idx])[::-1][: self.n]]] = True
        self.columns_ = cols[keep]
        return self

    def transform(self, X):
        return X[self.columns_]


class ColumnSelect(BaseEstimator, TransformerMixin):
    """Keep a fixed list of columns of a DataFrame."""

    def __init__(self, columns=None):
        self.columns = columns

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return X[list(self.columns)]


def _select(cols):
    return ColumnSelect(list(cols))


@dataclass
class ModelSpec:
    name: str
    columns: list[str]
    kind: str                      # "l1", "l2", "rf", "gbm"
    baseline: bool = False
    tune_C: bool = False


def make_pipeline(spec: ModelSpec, ds: Dataset, p3: dict, C: float = 1.0, seed: int = 0) -> Pipeline:
    wanted = set(spec.columns)
    blocks = {k: [c for c in v if c in wanted] for k, v in ds.blocks.items()}
    steps = [("select", _select(spec.columns)),
             ("topvar", BlockTopVariance(blocks, p3["features"]["top_genes_per_block"])),
             ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
             ("scale", StandardScaler())]
    m = p3["models"]
    if spec.kind == "l1":
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=C, class_weight="balanced", max_iter=2000)
    elif spec.kind == "l2":
        clf = LogisticRegression(penalty="l2", C=C, class_weight="balanced", max_iter=5000)
    elif spec.kind == "rf":
        clf = RandomForestClassifier(n_estimators=m["rf_trees"], max_features="sqrt", min_samples_leaf=2,
                                     class_weight="balanced", n_jobs=-1, random_state=seed)
    elif spec.kind == "gbm":
        clf = HistGradientBoostingClassifier(random_state=seed, class_weight="balanced", **m["gbm"])
    else:
        raise ValueError(spec.kind)
    steps.append(("clf", clf))
    return Pipeline(steps)


def model_specs(ds: Dataset) -> list[ModelSpec]:
    b = ds.blocks
    expr = [c for k, v in b.items() if k.endswith("expression") for c in v]
    props = b["cell-type proportions"]
    full = props + expr
    return [
        ModelSpec("Covariates only (batch + ancestry)", b["covariates"], "l2", baseline=True),
        ModelSpec("Cell-type proportions only", props, "l2", baseline=True),
        ModelSpec("Interferon score only", b["interferon scores"], "l2", baseline=True),
        ModelSpec("L1 logistic (all features)", full, "l1", tune_C=True),
        ModelSpec("Random forest (all features)", full, "rf"),
        ModelSpec("Gradient boosting (all features)", full, "gbm"),
    ]


# --------------------------------------------------------------------------
# cross-validation
# --------------------------------------------------------------------------
def _auc(y, p):
    return roc_auc_score(y, p) if len(np.unique(y)) == 2 else np.nan


def tune_C(spec, ds, tr, groups, p3, seed) -> float:
    grid = p3["cv"]["lr_C_grid"]
    inner = StratifiedGroupKFold(p3["cv"]["inner_folds"], shuffle=True, random_state=seed)
    Xtr, ytr, gtr = ds.X.iloc[tr], ds.y[tr], groups[tr]
    scores = []
    for C in grid:
        s = []
        for a, b in inner.split(Xtr, ytr, gtr):
            m = make_pipeline(spec, ds, p3, C=C, seed=seed).fit(Xtr.iloc[a], ytr[a])
            s.append(_auc(ytr[b], m.predict_proba(Xtr.iloc[b])[:, 1]))
        scores.append(np.nanmean(s))
    return float(grid[int(np.nanargmax(scores))])


def permute_block_auc_drop(model, Xte, yte, blocks, rng) -> dict[str, float]:
    base = _auc(yte, model.predict_proba(Xte)[:, 1])
    out = {}
    for name, cols in blocks.items():
        cols = [c for c in cols if c in Xte.columns]
        if not cols:
            continue
        Xp = Xte.copy()
        Xp[cols] = Xte[cols].to_numpy()[rng.permutation(len(Xte))]
        out[name] = base - _auc(yte, model.predict_proba(Xp)[:, 1])
    return out


@dataclass
class CVResult:
    preds: np.ndarray            # repeats x units
    chosen_C: list = field(default_factory=list)
    importance: list = field(default_factory=list)   # dicts per fold
    coefs: list = field(default_factory=list)        # pd.Series per fold (L1)


def run_cv(spec, ds, p3, seed=0, repeats=None, y=None, C_fixed=None, importance=False, verbose=True) -> CVResult:
    y = ds.y if y is None else y
    groups = ds.meta["donor"].to_numpy()
    reps = p3["cv"]["repeats"] if repeats is None else repeats
    preds = np.full((reps, len(y)), np.nan)
    res = CVResult(preds=preds)
    rng = np.random.default_rng(seed)
    imp_blocks = {k: v for k, v in ds.blocks.items() if k not in ("covariates", "interferon scores")}
    for r in range(reps):
        outer = StratifiedGroupKFold(p3["cv"]["outer_folds"], shuffle=True, random_state=seed + r)
        for k, (tr, te) in enumerate(outer.split(ds.X, y, groups)):
            C = 1.0
            if spec.tune_C:
                C = C_fixed if C_fixed is not None else tune_C(spec, _with_y(ds, y), tr, groups, p3, seed + 100 * r + k)
                res.chosen_C.append(C)
            m = make_pipeline(spec, ds, p3, C=C, seed=seed + r).fit(ds.X.iloc[tr], y[tr])
            preds[r, te] = m.predict_proba(ds.X.iloc[te])[:, 1]
            if importance:
                res.importance.append(permute_block_auc_drop(m, ds.X.iloc[te], y[te], imp_blocks, rng))
                clf = m.named_steps["clf"]
                if hasattr(clf, "coef_"):
                    res.coefs.append(pd.Series(clf.coef_.ravel(), index=m.named_steps["topvar"].columns_))
        if verbose:
            print(f"  {spec.name}: repeat {r + 1}/{reps}, pooled AUROC {_auc(y, preds[r]):.3f}", flush=True)
    return res


def _with_y(ds, y):
    return ds if y is ds.y else Dataset(X=ds.X, meta=ds.meta, y=y, blocks=ds.blocks)


# --------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------
def donor_level(meta: pd.DataFrame, y: np.ndarray, p: np.ndarray) -> pd.DataFrame:
    d = pd.DataFrame({"donor": meta["donor"].to_numpy(), "y": y, "p": p,
                      "ancestry": meta["ancestry"].to_numpy()})
    return d.groupby("donor").agg(y=("y", "first"), p=("p", "mean"), ancestry=("ancestry", "first"))


def bootstrap_auc(y, p, n, seed) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    vals = []
    for _ in range(n):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        vals.append(roc_auc_score(y[idx], p[idx]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def ifn_status(ds: Dataset, z_cut: float, ref_ct: str = "classical monocyte") -> pd.Series:
    """'IFN-high' / 'IFN-low' for SLE units, 'healthy' otherwise. z is relative
    to healthy units of the same batch, so batch shifts in the score cancel."""
    col = f"ifn|{ref_ct}"
    s = ds.X[col]
    out = pd.Series("healthy", index=ds.X.index)
    z = pd.Series(np.nan, index=ds.X.index)
    for b, idx in ds.meta.groupby("batch").groups.items():
        h = s.loc[idx][ds.meta.loc[idx, "disease_group"] == "healthy"].dropna()
        if len(h) >= 3 and h.std() > 0:
            z.loc[idx] = (s.loc[idx] - h.mean()) / h.std()
    sle = ds.meta["disease_group"] == "SLE"
    out[sle] = np.where(z[sle] > z_cut, "SLE IFN-high", "SLE IFN-low")
    out[sle & z.isna()] = "SLE IFN unknown"
    return out
