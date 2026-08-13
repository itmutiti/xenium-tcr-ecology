"""Unit tests for xenium_tcr_ecology.clone_ecology.apc_support (`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from xenium_tcr_ecology.clone_ecology.apc_support import (
    compute_neighbourhood_mean_score,
    compute_score_excess,
)


class TestComputeScoreExcess:
    def test_positive_excess_above_baseline(self):
        assert compute_score_excess(clone_value=0.3, section_baseline=-0.02) == pytest.approx(0.32)

    def test_negative_excess_below_baseline(self):
        assert compute_score_excess(clone_value=-0.3, section_baseline=-0.02) == pytest.approx(
            -0.28
        )

    def test_handles_negative_baseline_without_sign_flip(self):
        # A naive ratio (clone/baseline) with a negative denominator would
        # flip sign unpredictably; the difference must not.
        result = compute_score_excess(clone_value=0.1, section_baseline=-0.1)
        assert result == pytest.approx(0.2)

    def test_nan_baseline_gives_nan(self):
        assert np.isnan(compute_score_excess(clone_value=0.1, section_baseline=float("nan")))


class TestComputeNeighbourhoodMeanScore:
    def test_mean_of_real_neighbour_scores(self):
        # Cell 0 has two neighbours (score 2.0 and 4.0) -> mean 3.0.
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 0, 1, 2], [1, 2, 0, 0])), shape=(3, 3)
        )
        score = pd.Series([10.0, 2.0, 4.0], index=["c0", "c1", "c2"])
        result = compute_neighbourhood_mean_score(graph, score)
        assert result.loc["c0"] == 3.0

    def test_own_score_excluded_from_own_neighbourhood_mean(self):
        # Cell 0's own score (100.0) must not appear in its own result --
        # only its real neighbours' scores (no self-loops in the graph).
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        score = pd.Series([100.0, 5.0], index=["c0", "c1"])
        result = compute_neighbourhood_mean_score(graph, score)
        assert result.loc["c0"] == 5.0

    def test_zero_degree_cell_gets_nan(self):
        graph = sparse.csr_matrix((2, 2))
        score = pd.Series([1.0, 2.0], index=["c0", "c1"])
        result = compute_neighbourhood_mean_score(graph, score)
        assert result.loc["c0"] != result.loc["c0"]  # NaN

    def test_preserves_index(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        score = pd.Series([1.0, 2.0], index=["c0", "c1"])
        result = compute_neighbourhood_mean_score(graph, score)
        assert list(result.index) == ["c0", "c1"]
