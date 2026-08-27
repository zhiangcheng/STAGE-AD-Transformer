import numpy as np

from sexreg_ad.metrics import binary_metrics, expected_calibration_error


def test_probability_metrics_are_finite():
    labels = np.array([0, 0, 1, 1])
    probabilities = np.array([0.1, 0.3, 0.7, 0.9])
    metrics = binary_metrics(labels, probabilities)

    assert metrics["ScoreIsProbability"] is True
    assert metrics["ThresholdMode"] == "probability_0.5"
    assert metrics["AUROC"] == 1.0
    assert np.isfinite(metrics["Brier"])
    assert np.isfinite(expected_calibration_error(labels, probabilities))


def test_ranking_score_skips_probability_calibration():
    labels = np.array([0, 0, 1, 1])
    scores = np.array([2.0, 3.0, 7.0, 9.0])
    metrics = binary_metrics(labels, scores)

    assert metrics["ScoreIsProbability"] is False
    assert metrics["ThresholdMode"] == "top_prevalence_for_ranking_score"
    assert np.isnan(metrics["Brier"])
