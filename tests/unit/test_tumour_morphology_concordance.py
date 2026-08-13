"""Unit tests for xenium_tcr_ecology.tumour.morphology_concordance (`07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.tumour.morphology_concordance import (
    compute_malignancy_spatial_autocorrelation,
    select_extreme_regions,
)


class TestComputeMalignancySpatialAutocorrelation:
    def test_too_few_cells_returns_none(self):
        result = compute_malignancy_spatial_autocorrelation(
            np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([0.1, 0.9]), n_neighs=6
        )
        assert result["morans_i"] is None

    def test_spatially_clustered_pattern_has_higher_morans_i_than_random(self):
        rng = np.random.default_rng(0)
        n_per_cluster = 100

        # Two well-separated spatial clusters, each internally uniform in
        # malignancy_probability -- strong real spatial structure.
        x_clustered = np.concatenate(
            [rng.normal(0, 5, n_per_cluster), rng.normal(1000, 5, n_per_cluster)]
        )
        y_clustered = np.concatenate(
            [rng.normal(0, 5, n_per_cluster), rng.normal(0, 5, n_per_cluster)]
        )
        malignancy_clustered = np.concatenate(
            [np.full(n_per_cluster, 0.9), np.full(n_per_cluster, 0.1)]
        )

        # Same spatial layout, but malignancy values shuffled -- no real
        # spatial structure left.
        malignancy_random = rng.permutation(malignancy_clustered)

        clustered_result = compute_malignancy_spatial_autocorrelation(
            x_clustered, y_clustered, malignancy_clustered
        )
        random_result = compute_malignancy_spatial_autocorrelation(
            x_clustered, y_clustered, malignancy_random
        )

        assert clustered_result["morans_i"] > random_result["morans_i"]
        assert clustered_result["morans_i"] > 0.5


class TestSelectExtremeRegions:
    def test_selects_high_and_low_per_section(self):
        df = pd.DataFrame(
            {
                "section_id": ["S1"] * 10 + ["S2"] * 10,
                "malignancy_probability": list(np.linspace(0, 1, 10)) * 2,
            }
        )
        result = select_extreme_regions(df, n_per_side=2)
        assert set(result["section_id"].unique()) == {"S1", "S2"}
        s1 = result[result["section_id"] == "S1"]
        assert set(s1["region_type"]) == {"high", "low"}
        assert (
            s1[s1["region_type"] == "high"]["malignancy_probability"].min()
            > s1[s1["region_type"] == "low"]["malignancy_probability"].max()
        )
