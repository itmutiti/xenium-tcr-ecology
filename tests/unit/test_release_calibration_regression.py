"""Unit tests for xenium_tcr_ecology.release.calibration_regression (`17_statistical_closure_and_release/09_run_null_model_calibration_regression.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.release.calibration_regression import ci_overlaps


class TestCiOverlaps:
    def test_real_overlapping_intervals_return_true(self):
        assert ci_overlaps((0.0, 0.3), (0.1, 0.5)) is True

    def test_real_non_overlapping_intervals_return_false(self):
        assert ci_overlaps((0.0, 0.1), (0.5, 0.9)) is False

    def test_real_touching_intervals_return_true(self):
        assert ci_overlaps((0.0, 0.2), (0.2, 0.5)) is True
