"""Unit tests for xenium_tcr_ecology.graphs.synthetic_patterns (`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.graphs.synthetic_patterns import (
    assign_tumour_region,
    compute_domain_size_um,
    generate_background_points,
    sample_clone_cells,
)


class TestComputeDomainSizeUm:
    def test_returns_positive_finite_size(self):
        size = compute_domain_size_um(5000)
        assert size > 0
        assert np.isfinite(size)

    def test_more_cells_needs_a_larger_domain(self):
        assert compute_domain_size_um(10000) > compute_domain_size_um(5000)


class TestGenerateBackgroundPoints:
    def test_returns_requested_number_of_points_within_domain(self):
        rng = np.random.default_rng(0)
        points = generate_background_points(100, 500.0, rng)
        assert points.shape == (100, 2)
        assert (points >= 0).all() and (points <= 500.0).all()


class TestAssignTumourRegion:
    def test_hits_the_requested_fraction_approximately(self):
        rng = np.random.default_rng(0)
        points = generate_background_points(1000, 500.0, rng)
        mask = assign_tumour_region(points, tumour_fraction=0.1, rng=rng)
        assert mask.sum() == 100

    def test_tumour_region_is_spatially_contiguous(self):
        # A contiguous region's points should have a much smaller mean
        # pairwise distance than a random subset of the same size.
        rng = np.random.default_rng(0)
        points = generate_background_points(1000, 500.0, rng)
        mask = assign_tumour_region(points, tumour_fraction=0.1, rng=rng)
        tumour_points = points[mask]
        centre = tumour_points.mean(axis=0)
        max_dist_from_centre = np.linalg.norm(tumour_points - centre, axis=1).max()
        random_points = points[rng.choice(len(points), size=mask.sum(), replace=False)]
        random_max_dist = np.linalg.norm(random_points - random_points.mean(axis=0), axis=1).max()
        assert max_dist_from_centre < random_max_dist


class TestSampleCloneCells:
    def test_returns_requested_number_of_clone_cells(self):
        rng = np.random.default_rng(0)
        points = generate_background_points(500, 500.0, rng)
        is_tumour = assign_tumour_region(points, 0.1, rng)
        is_clone = sample_clone_cells(
            points, is_tumour, n_clone_cells=20, effect_size=1.0, length_scale=30.0, rng=rng
        )
        assert is_clone.sum() == 20

    def test_clone_cells_never_overlap_tumour_cells(self):
        rng = np.random.default_rng(0)
        points = generate_background_points(500, 500.0, rng)
        is_tumour = assign_tumour_region(points, 0.1, rng)
        is_clone = sample_clone_cells(
            points, is_tumour, n_clone_cells=20, effect_size=2.0, length_scale=30.0, rng=rng
        )
        assert not (is_clone & is_tumour).any()

    def test_positive_effect_size_places_clone_cells_closer_to_tumour_than_null(self):
        # The real, load-bearing property this module exists to guarantee:
        # a positive effect_size must produce clone cells that are, on
        # average, closer to the tumour region than a true null
        # (effect_size=0.0) would -- otherwise the "known, controllable
        # effect size" this milestone promises would not actually hold.
        from scipy.spatial import cKDTree

        rng_null = np.random.default_rng(42)
        points = generate_background_points(2000, 800.0, rng_null)
        is_tumour = assign_tumour_region(points, 0.073, rng_null)
        tumour_tree = cKDTree(points[is_tumour])

        rng_a = np.random.default_rng(1)
        null_clone = sample_clone_cells(
            points, is_tumour, 30, effect_size=0.0, length_scale=30.0, rng=rng_a
        )
        null_dist = tumour_tree.query(points[null_clone])[0].mean()

        rng_b = np.random.default_rng(2)
        effect_clone = sample_clone_cells(
            points, is_tumour, 30, effect_size=4.0, length_scale=30.0, rng=rng_b
        )
        effect_dist = tumour_tree.query(points[effect_clone])[0].mean()

        assert effect_dist < null_dist
