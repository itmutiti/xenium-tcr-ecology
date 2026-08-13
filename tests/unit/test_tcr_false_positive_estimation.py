"""Unit tests for xenium_tcr_ecology.tcr.false_positive_estimation (`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.tcr.false_positive_estimation import (
    compute_detection_spatial_autocorrelation,
    compute_non_tcell_detection_rate,
)


class TestComputeNonTcellDetectionRate:
    def test_computes_rate_for_intended_patient_only(self):
        counts = np.array([1, 0, 5, 0])
        patients = pd.Series(["P1", "P1", "P2", "P2"])
        result = compute_non_tcell_detection_rate(counts, patients, "P1")
        assert result["n_non_tcells"] == 2
        assert result["n_detected"] == 1
        assert result["detection_rate"] == 0.5

    def test_zero_cells_for_patient_returns_none_rate(self):
        counts = np.array([1, 1])
        patients = pd.Series(["P2", "P2"])
        result = compute_non_tcell_detection_rate(counts, patients, "P1")
        assert result["n_non_tcells"] == 0
        assert result["detection_rate"] is None


class TestComputeDetectionSpatialAutocorrelation:
    def test_too_few_cells_returns_none(self):
        result = compute_detection_spatial_autocorrelation(
            np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([1, 0]), n_neighs=6
        )
        assert result["morans_i"] is None

    def test_all_detected_or_none_detected_returns_none(self):
        x = np.arange(20, dtype=float)
        y = np.zeros(20)
        all_detected = np.ones(20, dtype=int)
        result = compute_detection_spatial_autocorrelation(x, y, all_detected, n_neighs=6)
        assert result["morans_i"] is None

        none_detected = np.zeros(20, dtype=int)
        result2 = compute_detection_spatial_autocorrelation(x, y, none_detected, n_neighs=6)
        assert result2["morans_i"] is None

    def test_spatially_clustered_detection_has_positive_morans_i(self):
        rng = np.random.default_rng(0)
        n_per_cluster = 30
        x = np.concatenate([rng.normal(0, 5, n_per_cluster), rng.normal(1000, 5, n_per_cluster)])
        y = np.concatenate([rng.normal(0, 5, n_per_cluster), rng.normal(0, 5, n_per_cluster)])
        detected = np.concatenate(
            [np.ones(n_per_cluster, dtype=int), np.zeros(n_per_cluster, dtype=int)]
        )
        result = compute_detection_spatial_autocorrelation(x, y, detected, n_neighs=6)
        assert result["morans_i"] > 0.5
