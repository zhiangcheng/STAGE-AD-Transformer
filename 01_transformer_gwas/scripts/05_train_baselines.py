#!/usr/bin/env python
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier

sys.path.append(str(Path(__file__).resolve().parents[1]))

from sexreg_ad.constants import VARIANT_FEATURES
from sexreg_ad.metrics import binary_metrics
from sexreg_ad.utils import ensure_dir


def parse_args():
    ap = argparse.ArgumentParser(
        description="Step 05: train and evaluate fixed-score and learned baselines."
    )
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--outdir", default="results/metrics")
    ap.add_argument("--max-train", type=int, default=500000)
    ap.add_argument("--seed", type=int, default=20260618)
    return ap.parse_args()


def main():
    args = parse_args()
    ensure_dir(args.outdir)

    train = pd.read_parquet(args.train)
    if len(train) > args.max_train:
        train = train.sample(args.max_train, random_state=args.seed)

    test = pd.read_parquet(args.test)

    cols = [
        c for c in test.columns
        if c in VARIANT_FEATURES
        or (c.startswith("gene") and not c.endswith("_name"))
        or c.startswith("tissue")
        or c.startswith("cell")
    ]

    Xtr = train[cols].replace([np.inf, -np.inf], np.nan).fillna(0)
    Xte = test[cols].replace([np.inf, -np.inf], np.nan).fillna(0)

    ytr = train.is_positive.astype(int).values
    yte = test.is_positive.astype(int).values

    rows = []

    score_defs = {
        "GWAS P-value ranking": test[
            ["neglog10P_total", "neglog10P_female", "neglog10P_male", "neglog10P_interaction"]
        ].max(axis=1).astype(float).values,
        "PIP-only": test[
            ["PIP_total", "PIP_female", "PIP_male", "PIP_interaction"]
        ].max(axis=1).astype(float).values,
        "Colocalization-only": test[
            ["coloc_brain_pph4", "coloc_blood_pph4", "coloc_immune_pph4"]
        ].max(axis=1).astype(float).values,
        "CADD/DeepSEA/Enformer-only": test[
            ["cadd_score", "deepsea_delta", "enformer_brain_delta", "sei_regulatory_score"]
        ].rank(pct=True).mean(axis=1).astype(float).values,
    }

    for name, score in score_defs.items():
        m = binary_metrics(yte, score)
        m["model"] = name
        rows.append(m)

    models = {
        "Elastic Net": make_pipeline(
            StandardScaler(),
            LogisticRegression(
                penalty="elasticnet",
                l1_ratio=0.15,
                solver="saga",
                class_weight="balanced",
                max_iter=1000,
                n_jobs=-1,
                random_state=args.seed,
            ),
        ),
        "Random Forest": RandomForestClassifier(
            n_estimators=200,
            max_depth=14,
            min_samples_leaf=20,
            class_weight="balanced_subsample",
            n_jobs=-1,
            random_state=args.seed,
        ),
        "MLP": make_pipeline(
            StandardScaler(),
            MLPClassifier(
                hidden_layer_sizes=(160, 80),
                max_iter=60,
                early_stopping=True,
                random_state=args.seed,
            ),
        ),
    }

    try:
        from xgboost import XGBClassifier
        spw = max(1, (len(ytr) - ytr.sum()) / max(1, ytr.sum()))
        models["XGBoost"] = XGBClassifier(
            n_estimators=400,
            max_depth=4,
            learning_rate=0.03,
            subsample=0.9,
            colsample_bytree=0.8,
            eval_metric="aucpr",
            n_jobs=-1,
            random_state=args.seed,
            scale_pos_weight=spw,
        )
    except Exception as e:
        print("Skip XGBoost", e)

    try:
        from lightgbm import LGBMClassifier
        models["LightGBM"] = LGBMClassifier(
            n_estimators=500,
            learning_rate=0.03,
            num_leaves=31,
            subsample=0.9,
            colsample_bytree=0.8,
            class_weight="balanced",
            n_jobs=-1,
            random_state=args.seed,
        )
    except Exception as e:
        print("Skip LightGBM", e)

    for name, model in models.items():
        print("Training", name)
        model.fit(Xtr, ytr)
        score = model.predict_proba(Xte)[:, 1]
        m = binary_metrics(yte, score)
        m["model"] = name
        rows.append(m)

    out = pd.DataFrame(rows).sort_values("AUPRC", ascending=False)
    out.to_csv(Path(args.outdir) / "baseline_metrics.csv", index=False)
    print(out)


if __name__ == "__main__":
    main()
