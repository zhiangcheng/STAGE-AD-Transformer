import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    f1_score,
    balanced_accuracy_score,
    brier_score_loss,
)
from sklearn.linear_model import LogisticRegression

from .constants import CLASS_NAMES
from .utils import topk_accuracy, mean_reciprocal_rank


def _as_clean_arrays(y_true, y_score):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    mask = np.isfinite(y_true) & np.isfinite(y_score)
    return y_true[mask], y_score[mask]


def _is_probability_score(y_score, tol=1e-8):
    y_score = np.asarray(y_score, dtype=float)
    if len(y_score) == 0:
        return False
    return np.nanmin(y_score) >= -tol and np.nanmax(y_score) <= 1.0 + tol


def _threshold_for_ranking_score(y_true, y_score):
    prevalence = float(np.mean(y_true))
    prevalence = min(max(prevalence, 1.0 / max(len(y_true), 1)), 0.999999)
    return np.quantile(y_score, 1.0 - prevalence)


def expected_calibration_error(y_true, y_prob, n_bins=10):
    y_true, y_prob = _as_clean_arrays(y_true, y_prob)
    if len(y_true) == 0 or not _is_probability_score(y_prob):
        return np.nan
    y_prob = np.clip(y_prob, 0.0, 1.0)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        lo, hi = bins[i], bins[i + 1]
        mask = (y_prob >= lo) & ((y_prob < hi) if i < n_bins - 1 else (y_prob <= hi))
        if mask.sum() == 0:
            continue
        ece += float(mask.mean()) * abs(float(y_true[mask].mean()) - float(y_prob[mask].mean()))
    return float(ece)


def calibration_slope(y_true, y_prob):
    y_true, y_prob = _as_clean_arrays(y_true, y_prob)
    if len(y_true) == 0 or len(np.unique(y_true)) < 2 or not _is_probability_score(y_prob):
        return np.nan
    y_prob = np.clip(y_prob, 1e-6, 1.0 - 1e-6)
    logit = np.log(y_prob / (1.0 - y_prob)).reshape(-1, 1)
    try:
        lr = LogisticRegression(solver="lbfgs", penalty=None, max_iter=1000)
        lr.fit(logit, y_true)
        return float(lr.coef_[0, 0])
    except TypeError:
        try:
            lr = LogisticRegression(solver="lbfgs", penalty="none", max_iter=1000)
            lr.fit(logit, y_true)
            return float(lr.coef_[0, 0])
        except Exception:
            return np.nan
    except Exception:
        return np.nan


def binary_metrics(y_true, y_score):
    """
    Robust binary metrics.

    AUROC/AUPRC accept any ranking score.
    Brier/ECE/calibration slope are computed only for scores in [0, 1].
    This prevents crashes when ranking-only baselines use -log10(P), CADD, etc.
    """
    y_true, y_score = _as_clean_arrays(y_true, y_score)
    out = {
        "AUROC": np.nan,
        "AUPRC": np.nan,
        "F1": np.nan,
        "BalancedAccuracy": np.nan,
        "Brier": np.nan,
        "ECE": np.nan,
        "CalibrationSlope": np.nan,
        "ScoreIsProbability": False,
        "ThresholdMode": "not_evaluated",
    }

    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return out

    try:
        out["AUROC"] = float(roc_auc_score(y_true, y_score))
    except Exception:
        pass

    try:
        out["AUPRC"] = float(average_precision_score(y_true, y_score))
    except Exception:
        pass

    is_prob = _is_probability_score(y_score)
    out["ScoreIsProbability"] = bool(is_prob)

    if is_prob:
        y_prob = np.clip(y_score, 0.0, 1.0)
        y_pred = (y_prob >= 0.5).astype(int)
        out["ThresholdMode"] = "probability_0.5"
        try:
            out["Brier"] = float(brier_score_loss(y_true, y_prob))
        except Exception:
            out["Brier"] = np.nan
        out["ECE"] = expected_calibration_error(y_true, y_prob)
        out["CalibrationSlope"] = calibration_slope(y_true, y_prob)
    else:
        threshold = _threshold_for_ranking_score(y_true, y_score)
        y_pred = (y_score >= threshold).astype(int)
        out["ThresholdMode"] = "top_prevalence_for_ranking_score"

    try:
        out["F1"] = float(f1_score(y_true, y_pred))
    except Exception:
        out["F1"] = np.nan

    try:
        out["BalancedAccuracy"] = float(balanced_accuracy_score(y_true, y_pred))
    except Exception:
        out["BalancedAccuracy"] = np.nan

    return out


def multiclass_metrics(y_true, prob):
    y_true = np.asarray(y_true).astype(int)
    prob = np.asarray(prob, dtype=float)
    pred = prob.argmax(axis=1)
    out = {
        "MacroF1": float(f1_score(y_true, pred, average="macro")),
        "WeightedF1": float(f1_score(y_true, pred, average="weighted")),
        "MulticlassBalancedAccuracy": float(balanced_accuracy_score(y_true, pred)),
    }
    for i, cls in enumerate(CLASS_NAMES):
        y_bin = (y_true == i).astype(int)
        try:
            out[f"AUPRC_{cls}"] = float(average_precision_score(y_bin, prob[:, i]))
        except Exception:
            out[f"AUPRC_{cls}"] = np.nan
    return out


def ranking_metrics(gene_y, gene_scores, tissue_y, tissue_scores, cell_y, cell_scores):
    return {
        "GeneTop1Acc": topk_accuracy(gene_y, gene_scores, k=1),
        "GeneTop3Acc": topk_accuracy(gene_y, gene_scores, k=3),
        "GeneMRR": mean_reciprocal_rank(gene_y, gene_scores),
        "TissueTop1Acc": topk_accuracy(tissue_y, tissue_scores, k=1),
        "TissueTop3Acc": topk_accuracy(tissue_y, tissue_scores, k=3),
        "CellTop1Acc": topk_accuracy(cell_y, cell_scores, k=1),
        "CellTop3Acc": topk_accuracy(cell_y, cell_scores, k=3),
    }


def topk_enrichment(df, score_col, y_col="is_positive", fracs=(0.001, 0.005, 0.01, 0.05)):
    rows = []
    global_rate = df[y_col].mean()
    for frac in fracs:
        k = max(1, int(len(df) * frac))
        top = df.sort_values(score_col, ascending=False).head(k)
        top_rate = top[y_col].mean()
        rows.append({
            "top_fraction": frac,
            "top_k": k,
            "top_positive_rate": top_rate,
            "global_positive_rate": global_rate,
            "fold_enrichment": top_rate / global_rate if global_rate > 0 else np.nan,
        })
    return pd.DataFrame(rows)
