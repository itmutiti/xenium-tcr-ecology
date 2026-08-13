"""Unit tests for xenium_tcr_ecology.validation.framework_generalisation (`16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.validation.framework_generalisation import subsample_real_points


class TestSubsampleRealPoints:
    def test_real_subsample_has_requested_size(self):
        rng = np.random.default_rng(0)
        points = np.column_stack([np.arange(1000.0), np.zeros(1000)])
        result = subsample_real_points(points, n_cells=100, rng=rng)
        assert len(result) == 100

    def test_real_subsample_points_are_unique(self):
        rng = np.random.default_rng(1)
        points = np.column_stack([np.arange(1000.0), np.zeros(1000)])
        result = subsample_real_points(points, n_cells=100, rng=rng)
        assert len(np.unique(result, axis=0)) == 100

    def test_real_subsample_never_exceeds_real_population_size(self):
        rng = np.random.default_rng(2)
        points = np.column_stack([np.arange(50.0), np.zeros(50)])
        result = subsample_real_points(points, n_cells=200, rng=rng)
        assert len(result) == 50

    def test_real_subsample_is_spatially_contiguous_not_scattered(self):
        # Two real, well-separated clusters -- a real spatially-
        # contiguous subsample of one cluster's own size should stay
        # entirely within that one cluster, not draw from both.
        rng = np.random.default_rng(3)
        cluster_a = np.random.default_rng(10).normal(0, 1, size=(50, 2))
        cluster_b = np.random.default_rng(11).normal(1000, 1, size=(50, 2))
        points = np.vstack([cluster_a, cluster_b])
        result = subsample_real_points(points, n_cells=50, rng=rng)
        spans_x = result[:, 0].max() - result[:, 0].min()
        assert spans_x < 100  # real single-cluster span, not the real ~1000 cross-cluster span
