"""Unit tests for xenium_tcr_ecology.graphs.null_models (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.graphs.null_models import (
    clopper_pearson_ci,
    compute_mean_distance_to_tumour,
    degree_matched_permutation_pvalue,
    deterministic_seed,
    free_permutation_pvalue,
    graph_distance_permutation_pvalue,
)
from xenium_tcr_ecology.graphs.candidate_graphs import build_radius_graph
from xenium_tcr_ecology.graphs.synthetic_patterns import assign_tumour_region, sample_clone_cells


def _make_attracted_pattern(rng, n=300, attraction=True):
    # Uses the same real, contiguous tumour-region construction and
    # importance-weighted clone sampling as production code (`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`),
    # not an ad hoc test-only fixture -- a scattered (non-contiguous)
    # "tumour" subset with clone cells chosen by distance-to-CENTROID
    # would not actually be close to any individual tumour point, and
    # produced a real false failure here before this fix.
    points = rng.uniform(0, 500, size=(n, 2))
    is_tumour = assign_tumour_region(points, tumour_fraction=0.1, rng=rng)
    effect_size = 4.0 if attraction else 0.0
    is_clone = sample_clone_cells(
        points, is_tumour, n_clone_cells=20, effect_size=effect_size, length_scale=30.0, rng=rng
    )
    return points, is_tumour, is_clone


class TestComputeMeanDistanceToTumour:
    def test_zero_for_colocated_clone_and_tumour(self):
        points = np.array([[0.0, 0.0], [0.0, 0.0], [10.0, 10.0]])
        is_tumour = np.array([True, False, False])
        is_clone = np.array([False, True, False])
        assert compute_mean_distance_to_tumour(points, is_tumour, is_clone) == 0.0


class TestFreePermutationPvalue:
    def test_attracted_pattern_gives_low_pvalue(self):
        rng = np.random.default_rng(0)
        points, is_tumour, is_clone = _make_attracted_pattern(rng, attraction=True)
        p = free_permutation_pvalue(points, is_tumour, is_clone, rng, n_permutations=199)
        assert p < 0.05

    def test_random_pattern_gives_non_significant_pvalue_typically(self):
        rng = np.random.default_rng(1)
        points, is_tumour, is_clone = _make_attracted_pattern(rng, attraction=False)
        p = free_permutation_pvalue(points, is_tumour, is_clone, rng, n_permutations=199)
        assert p > 0.05

    def test_pvalue_is_in_valid_range(self):
        rng = np.random.default_rng(2)
        points, is_tumour, is_clone = _make_attracted_pattern(rng, attraction=True)
        p = free_permutation_pvalue(points, is_tumour, is_clone, rng, n_permutations=199)
        assert 0.0 < p <= 1.0


class TestDegreeMatchedPermutationPvalue:
    def test_attracted_pattern_gives_low_pvalue(self):
        rng = np.random.default_rng(0)
        points, is_tumour, is_clone = _make_attracted_pattern(rng, attraction=True)
        graph = build_radius_graph(points, radius_um=50.0)
        degree = np.asarray(graph.sum(axis=1)).ravel()
        import pandas as pd

        strata = pd.qcut(degree, q=4, labels=False, duplicates="drop")
        p = degree_matched_permutation_pvalue(
            points, is_tumour, is_clone, strata, rng, n_permutations=199
        )
        assert p < 0.1


class TestGraphDistancePermutationPvalue:
    def test_attracted_pattern_gives_low_pvalue(self):
        rng = np.random.default_rng(0)
        points, is_tumour, is_clone = _make_attracted_pattern(rng, attraction=True)
        graph = build_radius_graph(points, radius_um=50.0)
        p = graph_distance_permutation_pvalue(graph, is_tumour, is_clone, rng, n_permutations=199)
        assert p < 0.1


class TestClopperPearsonCi:
    def test_zero_successes_gives_lower_bound_zero(self):
        lower, upper = clopper_pearson_ci(0, 10)
        assert lower == 0.0
        assert upper > 0.0

    def test_all_successes_gives_upper_bound_one(self):
        lower, upper = clopper_pearson_ci(10, 10)
        assert upper == 1.0
        assert lower < 1.0

    def test_interval_contains_the_point_estimate(self):
        lower, upper = clopper_pearson_ci(1, 10)
        assert lower <= 0.1 <= upper


class TestDeterministicSeed:
    def test_real_same_inputs_give_same_seed_within_a_process(self):
        assert deterministic_seed("SYN00", 0.0, "calibration") == deterministic_seed(
            "SYN00", 0.0, "calibration"
        )

    def test_real_different_inputs_give_different_seeds(self):
        assert deterministic_seed("SYN00", 0.0, "calibration") != deterministic_seed(
            "SYN01", 0.0, "calibration"
        )
        assert deterministic_seed("SYN00", 0.0, "calibration") != deterministic_seed(
            "SYN00", 0.1, "calibration"
        )

    def test_real_seed_is_a_valid_numpy_default_rng_seed(self):
        seed = deterministic_seed("SYN00", 0.0, "calibration")
        np.random.default_rng(seed)  # should not raise
        assert 0 <= seed < 2**32

    def test_real_stable_across_separate_real_process_invocations(self):
        # A real regression test for the exact bug `17_statistical_closure_and_release/09_run_null_model_calibration_regression.py` found:
        # Python's built-in hash() randomises str-containing tuples per
        # process; deterministic_seed must not.
        import subprocess
        import sys

        code = "from xenium_tcr_ecology.graphs.null_models import deterministic_seed; print(deterministic_seed('SYN00', 0.0, 'calibration'))"
        first = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        ).stdout.strip()
        second = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        ).stdout.strip()
        assert first == second
