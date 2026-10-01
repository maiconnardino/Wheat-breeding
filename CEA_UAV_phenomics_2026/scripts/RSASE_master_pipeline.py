"""
RSASE UAV remote-sensing reproducibility workflow
Manuscript: Longitudinal UAV remote sensing for wheat breeding decisions:
            cross-experiment transferability and early selection utility
Target journal: Remote Sensing Applications: Society and Environment

This Python workflow is the production reference for the RSASE revision.
All data-dependent preprocessing is fitted within training data only.

Fixed longitudinal XGBoost specification (frozen for the RSASE revision):
    n_estimators=120
    max_depth=2
    learning_rate=0.05
    subsample=0.80
    colsample_bytree=0.80
    min_child_weight=1
    reg_lambda=1.0
    reg_alpha=0.0
    gamma=0.0
    objective="reg:squarederror"
    random_state=20260922
    n_jobs=1
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.compose import TransformedTargetRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GridSearchCV, GroupKFold, KFold, RepeatedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from xgboost import XGBRegressor

SEED = 20260922
DATA_FILE = Path("data/CEA_wheat_phenomics_154_genotypes.csv")
OUT = Path("results")
OUT.mkdir(exist_ok=True)

XGB_LONGITUDINAL_PARAMS = dict(
    n_estimators=120,
    max_depth=2,
    learning_rate=0.05,
    subsample=0.80,
    colsample_bytree=0.80,
    min_child_weight=1,
    reg_lambda=1.0,
    reg_alpha=0.0,
    gamma=0.0,
    objective="reg:squarederror",
    random_state=SEED,
    n_jobs=1,
    verbosity=0,
)

RF_GRID = [
    {"model__n_estimators": [100], "model__max_features": [0.3], "model__min_samples_leaf": [1]},
    {"model__n_estimators": [100], "model__max_features": [0.7], "model__min_samples_leaf": [1]},
    {"model__n_estimators": [100], "model__max_features": [1.0], "model__min_samples_leaf": [3]},
]
SVM_GRID = [
    {"model__C": [1.0], "model__epsilon": [0.1], "model__gamma": ["scale"]},
    {"model__C": [10.0], "model__epsilon": [0.1], "model__gamma": ["scale"]},
    {"model__C": [100.0], "model__epsilon": [0.1], "model__gamma": ["scale"]},
]
XGB_GRID = [
    {"model__n_estimators": [120], "model__max_depth": [2], "model__learning_rate": [0.05],
     "model__subsample": [0.8], "model__colsample_bytree": [0.8]},
    {"model__n_estimators": [120], "model__max_depth": [3], "model__learning_rate": [0.05],
     "model__subsample": [0.8], "model__colsample_bytree": [0.8]},
    {"model__n_estimators": [120], "model__max_depth": [2], "model__learning_rate": [0.10],
     "model__subsample": [0.8], "model__colsample_bytree": [1.0]},
]

RIDGE_ALPHAS = np.logspace(-3, 3, 25)
ENET_ALPHA = [0.001, 0.01, 0.1, 1.0, 10.0]
ENET_L1 = [0.1, 0.5, 0.9]


def read_data(path: Path = DATA_FILE) -> pd.DataFrame:
    d = pd.read_csv(path)
    assert len(d) == 154
    assert d["GEN"].nunique() == 154
    assert {"GEN", "Genetic material", "DH", "GY"}.issubset(d.columns)
    d["GEN_num"] = d["GEN"].str.replace("^G", "", regex=True).astype(int)
    d["Stage"] = np.select(
        [d.GEN_num <= 49, d.GEN_num <= 105, d.GEN_num <= 124],
        ["E1", "E2", "E3"], default="E4"
    )
    assert d.Stage.value_counts().to_dict() == {"E2": 56, "E1": 49, "E4": 30, "E3": 19}
    return d


def flight_columns(d: pd.DataFrame, max_flight: int = 10) -> list[str]:
    out = []
    for c in d.columns:
        m = re.search(r"\.flight(\d{2})$", c)
        if m and int(m.group(1)) <= max_flight:
            out.append(c)
    if max_flight == 10:
        assert len(out) == 120
    return out


class CorrFilter:
    """Train-only greedy absolute-correlation filter."""
    def __init__(self, threshold: float = 0.99):
        self.threshold = threshold
        self.keep_: list[str] | None = None

    def fit(self, X: pd.DataFrame):
        corr = X.corr().abs()
        upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
        drop = [c for c in upper.columns if (upper[c] > self.threshold).any()]
        self.keep_ = [c for c in X.columns if c not in drop]
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        if self.keep_ is None:
            raise RuntimeError("CorrFilter must be fitted first.")
        return X.loc[:, self.keep_]


def regression_metrics(y: np.ndarray, p: np.ndarray) -> dict:
    y = np.asarray(y, float); p = np.asarray(p, float)
    pearson = float(np.corrcoef(y, p)[0, 1])
    rho = float(spearmanr(y, p).statistic)
    r2 = float(1.0 - np.sum((y - p) ** 2) / np.sum((y - np.mean(y)) ** 2))
    return {
        "Pearson_r": pearson,
        "Spearman_rho": rho,
        "R2": r2,
        "RMSE": float(np.sqrt(mean_squared_error(y, p))),
        "MAE": float(mean_absolute_error(y, p)),
    }


def cross_experiment_from_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    """Add Spearman rho to already audited held-out-domain predictions."""
    rows = []
    for (stage, trait, model), g in predictions.groupby(["TargetStage", "Trait", "Model"]):
        m = regression_metrics(g["truth"].to_numpy(), g["pred"].to_numpy())
        rows.append({"TargetStage": stage, "Trait": trait, "Model": model, "n": len(g), **m})
    return pd.DataFrame(rows)


def fixed_longitudinal_xgb(
    d: pd.DataFrame,
    stage: str,
    trait: str,
    max_flight: int,
    split_membership: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    """
    Refit the frozen RSASE longitudinal XGBoost specification using the archived
    repeated-fourfold split membership. Correlation filtering is fitted inside
    each training fold. Predictions are averaged across repeated OOF appearances.
    """
    dat = d.loc[d.Stage.eq(stage)].copy()
    cols = flight_columns(d, max_flight)
    acc = {g: [] for g in dat.GEN}
    split_rows = []

    for split_id, sg in split_membership.groupby("split_id"):
        test_ids = set(sg.GEN)
        tr = dat.loc[~dat.GEN.isin(test_ids)].copy()
        te = dat.loc[dat.GEN.isin(test_ids)].copy()

        cf = CorrFilter(0.99).fit(tr[cols])
        keep = cf.keep_
        model = XGBRegressor(**XGB_LONGITUDINAL_PARAMS)
        model.fit(tr[keep], tr[trait])
        pred = model.predict(te[keep])

        for gen, truth, phat in zip(te.GEN, te[trait], pred):
            acc[gen].append(float(phat))
            split_rows.append({
                "Stage": stage, "Trait": trait, "flight_max": max_flight,
                "Model": "XGBoost", "split_id": int(split_id), "GEN": gen,
                "truth": float(truth), "pred": float(phat),
                "n_predictors_after_filter": len(keep),
            })

    oof = pd.DataFrame({
        "Stage": stage,
        "Trait": trait,
        "flight_max": max_flight,
        "Model": "XGBoost",
        "GEN": dat.GEN,
        "truth": dat[trait].astype(float),
        "pred": [float(np.mean(acc[g])) for g in dat.GEN],
    })
    return pd.DataFrame(split_rows), regression_metrics(oof.truth.to_numpy(), oof.pred.to_numpy())


def top_k_summary(truth: pd.Series, pred: pd.Series, q: float) -> dict:
    k = int(math.ceil(len(truth) * q))
    obs = set(truth.nlargest(k).index)
    est = set(pred.nlargest(k).index)
    recovery = len(obs & est) / k
    return {
        "selection_fraction": q,
        "k": k,
        "recovery": recovery,
        "random_expectation": k / len(truth),
        "enrichment": recovery / (k / len(truth)),
    }


if __name__ == "__main__":
    d = read_data()
    print("Data audit passed:", d.shape)
    print("Frozen longitudinal XGBoost parameters:")
    for k, v in XGB_LONGITUDINAL_PARAMS.items():
        print(f"  {k} = {v}")
    print(
        "\nCross-experiment Spearman rho is calculated from the final held-out-domain "
        "prediction files so that ranking transfer is reported separately from Pearson r "
        "and absolute calibration."
    )
