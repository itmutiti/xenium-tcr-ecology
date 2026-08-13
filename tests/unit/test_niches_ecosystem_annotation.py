"""Unit tests for xenium_tcr_ecology.niches.ecosystem_annotation (`10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.niches.ecosystem_annotation import (
    compute_enrichment_ratios,
    label_ecosystem,
)


class TestComputeEnrichmentRatios:
    def test_ratio_is_centroid_over_global_mean(self):
        centroids = pd.DataFrame({"T_cell": [0.4, 0.1], "Fibroblast": [0.1, 0.5]}, index=[1, 2])
        global_mean = pd.Series({"T_cell": 0.2, "Fibroblast": 0.1})
        result = compute_enrichment_ratios(centroids, global_mean, ["T_cell", "Fibroblast"])
        assert result.loc[1, "T_cell"] == 2.0
        assert result.loc[1, "Fibroblast"] == 1.0
        assert result.loc[2, "T_cell"] == 0.5
        assert result.loc[2, "Fibroblast"] == 5.0


class TestLabelEcosystem:
    def test_single_lineage_clears_threshold(self):
        row = pd.Series({"Fibroblast": 4.7, "T_cell": 0.8, "Endothelial": 0.6})
        assert label_ecosystem(row, threshold=2.0) == "Fibroblast niche"

    def test_two_lineages_clear_threshold_ordered_by_enrichment_descending(self):
        row = pd.Series({"Perivascular_SmoothMuscle": 5.5, "Endothelial": 4.4, "Fibroblast": 1.2})
        assert label_ecosystem(row, threshold=2.0) == "Perivascular/Endothelial niche"

    def test_no_lineage_clears_threshold_gives_mixed_label(self):
        row = pd.Series({"Mast_cell": 1.4, "Epithelial_Tumour": 1.25, "T_cell": 0.8})
        assert label_ecosystem(row, threshold=2.0) == "Mixed/non-specific niche"

    def test_three_way_tie_all_included(self):
        row = pd.Series({"NK_cell": 2.9, "T_cell": 2.7, "Myeloid": 2.6, "B_cell": 1.5})
        assert label_ecosystem(row, threshold=2.0) == "NK-cell/T-cell/Myeloid niche"
