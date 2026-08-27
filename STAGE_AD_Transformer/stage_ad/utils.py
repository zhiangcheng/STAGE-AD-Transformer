"""Small reusable numerical and filesystem helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.stats import norm


def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def p_from_z(z_score: np.ndarray) -> np.ndarray:
    return np.clip(2 * norm.sf(np.abs(z_score)), 1e-300, 1.0)


def neglog10(p_value: np.ndarray) -> np.ndarray:
    return -np.log10(np.clip(p_value, 1e-300, 1.0))


def chromosome_order() -> list[str]:
    return [str(index) for index in range(1, 23)] + ["X"]


def topk_accuracy(labels: np.ndarray, scores: np.ndarray, k: int = 1) -> float:
    labels = np.asarray(labels)
    order = np.argsort(-scores, axis=1)[:, :k]
    return float(np.mean([labels[index] in order[index] for index in range(len(labels))]))


def mean_reciprocal_rank(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels)
    order = np.argsort(-scores, axis=1)
    reciprocal_ranks = []
    for index in range(len(labels)):
        position = np.where(order[index] == labels[index])[0]
        reciprocal_ranks.append(1 / (position[0] + 1) if len(position) else 0)
    return float(np.mean(reciprocal_ranks))
