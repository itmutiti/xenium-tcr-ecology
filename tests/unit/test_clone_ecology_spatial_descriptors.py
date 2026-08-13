"""Unit tests for xenium_tcr_ecology.clone_ecology.spatial_descriptors (`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.clone_ecology.spatial_descriptors import (
    compute_border_enrichment,
    compute_clark_evans_index,
    compute_convex_hull_area,
    compute_dispersion,
    compute_rarefied_clark_evans,
    compute_rarefied_dispersion,
    compute_rarefied_domain_richness,
    filter_to_primary_cohort,
)


class TestFilterToPrimaryCohort:
    def test_excludes_non_primary_sections(self):
        manifest = pd.DataFrame(
            {"section_id": ["A", "B", "C"], "included_in_primary_hnscc_cohort": [True, False, True]}
        )
        result = filter_to_primary_cohort(manifest)
        assert list(result["section_id"]) == ["A", "C"]


class TestComputeDispersion:
    def test_zero_for_identical_points(self):
        x = np.array([5.0, 5.0, 5.0])
        y = np.array([3.0, 3.0, 3.0])
        assert compute_dispersion(x, y) == 0.0

    def test_positive_for_spread_points(self):
        x = np.array([0.0, 10.0])
        y = np.array([0.0, 0.0])
        assert compute_dispersion(x, y) == 5.0


class TestComputeConvexHullArea:
    def test_unit_square_area_is_one(self):
        x = np.array([0.0, 1.0, 1.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        assert np.isclose(compute_convex_hull_area(x, y), 1.0)

    def test_fewer_than_three_points_returns_none(self):
        assert compute_convex_hull_area(np.array([0.0, 1.0]), np.array([0.0, 1.0])) is None


class TestComputeClarkEvansIndex:
    def test_tight_cluster_gives_index_below_one(self):
        # Points much closer together than CSR expectation at this density -> clustered.
        rng = np.random.default_rng(0)
        x = rng.normal(0, 0.01, 20)
        y = rng.normal(0, 0.01, 20)
        result = compute_clark_evans_index(x, y, reference_density=1.0)
        assert result < 1.0

    def test_single_point_gives_nan(self):
        assert np.isnan(
            compute_clark_evans_index(np.array([0.0]), np.array([0.0]), reference_density=1.0)
        )

    def test_zero_density_gives_nan(self):
        assert np.isnan(
            compute_clark_evans_index(
                np.array([0.0, 1.0]), np.array([0.0, 1.0]), reference_density=0.0
            )
        )


class TestComputeBorderEnrichment:
    def test_ratio_of_two_over_baseline(self):
        # Section overall: 25% margin. Clone: 50% margin -> enrichment 2x.
        section = pd.Series(["inner_margin", "distal_stroma", "distal_stroma", "distal_stroma"])
        clone = pd.Series(["inner_margin", "distal_stroma"])
        result = compute_border_enrichment(clone, section)
        assert np.isclose(result, 2.0)

    def test_zero_baseline_gives_nan(self):
        section = pd.Series(["distal_stroma", "distal_stroma"])
        clone = pd.Series(["distal_stroma"])
        assert np.isnan(compute_border_enrichment(clone, section))


class TestRarefaction:
    def test_rarefied_dispersion_is_nan_below_n_rarefy(self):
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 1.0])
        rng = np.random.default_rng(1)
        result = compute_rarefied_dispersion(x, y, n_rarefy=5, n_iterations=10, rng=rng)
        assert np.isnan(result)

    def test_rarefied_dispersion_computed_when_enough_cells(self):
        rng = np.random.default_rng(1)
        x = rng.normal(0, 1, 20)
        y = rng.normal(0, 1, 20)
        result = compute_rarefied_dispersion(x, y, n_rarefy=5, n_iterations=50, rng=rng)
        assert result > 0
        assert not np.isnan(result)

    def test_rarefied_clark_evans_is_nan_below_n_rarefy(self):
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.0, 1.0, 2.0])
        rng = np.random.default_rng(1)
        result = compute_rarefied_clark_evans(
            x, y, reference_density=1.0, n_rarefy=5, n_iterations=10, rng=rng
        )
        assert np.isnan(result)

    def test_rarefied_domain_richness_counts_distinct_domains(self):
        # 10 cells, only 2 distinct domains -- rarefied richness must be <= 2.
        domain_ids = np.array([1] * 5 + [2] * 5)
        rng = np.random.default_rng(2)
        result = compute_rarefied_domain_richness(domain_ids, n_rarefy=5, n_iterations=100, rng=rng)
        assert 1.0 <= result <= 2.0

    def test_rarefied_domain_richness_below_n_rarefy_is_nan(self):
        domain_ids = np.array([1, 2, 3])
        rng = np.random.default_rng(2)
        result = compute_rarefied_domain_richness(domain_ids, n_rarefy=5, n_iterations=10, rng=rng)
        assert np.isnan(result)
