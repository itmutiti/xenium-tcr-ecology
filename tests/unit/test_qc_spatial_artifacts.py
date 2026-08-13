"""Unit tests for xenium_tcr_ecology.qc.spatial_artifacts (`04_quality_control/02_detect_spatial_qc_artifacts.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.qc.spatial_artifacts import MAD_FLAG_THRESHOLD, _modified_z_scores


class TestModifiedZScores:
    def test_flags_a_real_outlier(self):
        values = pd.Series([10.0, 10.5, 9.8, 10.2, 9.9, 50.0])  # last value is a clear outlier
        scores = _modified_z_scores(values)
        assert abs(scores.iloc[-1]) > MAD_FLAG_THRESHOLD
        assert all(abs(s) < MAD_FLAG_THRESHOLD for s in scores.iloc[:-1])

    def test_zero_mad_returns_zero_scores_not_nan_or_inf(self):
        """All-identical values would divide by zero MAD -- must not crash
        or produce NaN/inf, which would silently break the flagging logic."""
        values = pd.Series([5.0, 5.0, 5.0, 5.0])
        scores = _modified_z_scores(values)
        assert (scores == 0).all()

    def test_no_false_flags_on_tight_cluster(self):
        values = pd.Series([10.0, 10.1, 9.9, 10.05, 9.95, 10.02])
        scores = _modified_z_scores(values)
        assert all(abs(s) < MAD_FLAG_THRESHOLD for s in scores)
