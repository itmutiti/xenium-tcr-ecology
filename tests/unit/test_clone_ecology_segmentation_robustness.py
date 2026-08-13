"""Unit tests for xenium_tcr_ecology.clone_ecology.segmentation_robustness (`13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.clone_ecology.segmentation_robustness import compare_resegmented_vs_rest


class TestCompareResegmentedVsRest:
    def test_real_mean_split_by_section_membership(self):
        engagement = pd.DataFrame(
            {
                "section_id": ["P19_run1", "P19_run1", "P01_run1", "P01_run1"],
                "malignant_adjacency": [0.4, 0.6, 0.1, 0.3],
            }
        )
        result = compare_resegmented_vs_rest(
            engagement, resegmented_sections=["P19_run1"], metrics=["malignant_adjacency"]
        )
        row = result.iloc[0]
        assert row["mean_in_subset"] == 0.5
        assert row["mean_outside_subset"] == 0.2
        assert row["n_in_subset"] == 2
        assert row["n_outside_subset"] == 2

    def test_missing_section_gives_zero_count_not_error(self):
        engagement = pd.DataFrame({"section_id": ["P01_run1"], "malignant_adjacency": [0.5]})
        result = compare_resegmented_vs_rest(
            engagement, resegmented_sections=["P19_run1"], metrics=["malignant_adjacency"]
        )
        assert result.iloc[0]["n_in_subset"] == 0

    def test_multiple_metrics_each_get_own_row(self):
        engagement = pd.DataFrame(
            {
                "section_id": ["P19_run1", "P01_run1"],
                "malignant_adjacency": [0.5, 0.2],
                "engagement_ratio": [1.2, 0.8],
            }
        )
        result = compare_resegmented_vs_rest(
            engagement,
            resegmented_sections=["P19_run1"],
            metrics=["malignant_adjacency", "engagement_ratio"],
        )
        assert set(result["metric"]) == {"malignant_adjacency", "engagement_ratio"}
