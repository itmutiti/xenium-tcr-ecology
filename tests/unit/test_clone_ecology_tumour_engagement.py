"""Unit tests for xenium_tcr_ecology.clone_ecology.tumour_engagement (`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.clone_ecology.tumour_engagement import (
    compute_engagement_ratio,
    compute_penetration,
)


class TestComputeEngagementRatio:
    def test_above_opportunity_baseline_gives_ratio_above_one(self):
        result = compute_engagement_ratio(clone_mean_adjacency=0.4, section_mean_adjacency=0.2)
        assert result == 2.0

    def test_below_opportunity_baseline_gives_ratio_below_one(self):
        result = compute_engagement_ratio(clone_mean_adjacency=0.1, section_mean_adjacency=0.2)
        assert result == 0.5

    def test_zero_baseline_gives_nan(self):
        assert np.isnan(
            compute_engagement_ratio(clone_mean_adjacency=0.1, section_mean_adjacency=0.0)
        )

    def test_nan_baseline_gives_nan(self):
        assert np.isnan(
            compute_engagement_ratio(clone_mean_adjacency=0.1, section_mean_adjacency=float("nan"))
        )


class TestComputePenetration:
    def test_all_outside_tumour_gives_zero_fraction(self):
        distances = pd.Series([5.0, 10.0, 20.0])
        result = compute_penetration(distances)
        assert result["fraction_inside_tumour"] == 0.0
        assert result["min_signed_distance_um"] == 5.0

    def test_some_inside_tumour(self):
        distances = pd.Series([-3.0, 5.0, -1.0, 10.0])
        result = compute_penetration(distances)
        assert result["fraction_inside_tumour"] == 0.5
        assert result["min_signed_distance_um"] == -3.0

    def test_all_null_gives_nan(self):
        distances = pd.Series([np.nan, np.nan])
        result = compute_penetration(distances)
        assert np.isnan(result["fraction_inside_tumour"])
        assert np.isnan(result["min_signed_distance_um"])

    def test_null_values_excluded_from_computation(self):
        distances = pd.Series([-2.0, np.nan, 4.0])
        result = compute_penetration(distances)
        assert result["fraction_inside_tumour"] == 0.5
        assert result["min_signed_distance_um"] == -2.0
