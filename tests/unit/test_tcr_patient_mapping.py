"""Unit tests for xenium_tcr_ecology.tcr.patient_mapping (`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.tcr.patient_mapping import (
    compute_probe_patient_detection,
    evaluate_patient_specificity,
)


class TestComputeProbePatientDetection:
    def test_computes_detection_rate_per_candidate_patient(self):
        counts = np.array([1, 0, 0, 2, 0])
        patients = pd.Series(["P1", "P1", "P2", "P2", "P2"])
        result = compute_probe_patient_detection(counts, patients, ["P1", "P2"])
        assert result.loc["P1", "n_tcells"] == 2
        assert result.loc["P1", "n_detected"] == 1
        assert result.loc["P1", "detection_rate"] == 0.5
        assert result.loc["P2", "detection_rate"] == pytest.approx(1 / 3)

    def test_excludes_non_candidate_patients(self):
        counts = np.array([5, 5])
        patients = pd.Series(["P1", "P99"])
        result = compute_probe_patient_detection(counts, patients, ["P1"])
        assert list(result.index) == ["P1"]


class TestTestPatientSpecificity:
    def test_clear_single_patient_signal_is_significant(self):
        summary = pd.DataFrame(
            {
                "n_tcells": [1000, 1000, 1000],
                "n_detected": [500, 2, 3],
                "detection_rate": [0.5, 0.002, 0.003],
            },
            index=["P1", "P2", "P3"],
        )
        result = evaluate_patient_specificity(summary)
        assert result["top_patient"] == "P1"
        assert result["pvalue"] < 0.001
        assert result["direction_consistent"]

    def test_no_differentiation_gives_high_pvalue(self):
        summary = pd.DataFrame(
            {
                "n_tcells": [100, 100, 100],
                "n_detected": [5, 5, 5],
                "detection_rate": [0.05, 0.05, 0.05],
            },
            index=["P1", "P2", "P3"],
        )
        result = evaluate_patient_specificity(summary)
        assert result["pvalue"] > 0.05

    def test_single_candidate_patient_cannot_be_tested(self):
        summary = pd.DataFrame(
            {"n_tcells": [100], "n_detected": [5], "detection_rate": [0.05]}, index=["P1"]
        )
        result = evaluate_patient_specificity(summary)
        assert result["pvalue"] is None

    def test_zero_cells_for_all_candidates_cannot_be_tested(self):
        summary = pd.DataFrame(
            {"n_tcells": [0, 0], "n_detected": [0, 0], "detection_rate": [np.nan, np.nan]},
            index=["P1", "P2"],
        )
        result = evaluate_patient_specificity(summary)
        assert result["top_patient"] is None
