"""Unit tests for xenium_tcr_ecology.interactions.prioritisation (`14_spatial_interactions_and_barriers/05_prioritise_testable_interactions.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.interactions.prioritisation import (
    aggregate_spatial_scores,
    compute_composite_priority,
)


class TestAggregateSpatialScores:
    def test_real_spatial_specificity_is_fraction_of_significant_sections(self):
        scores = pd.DataFrame(
            {
                "section_id": ["S1", "S2", "S3", "S4"],
                "sender_receiver_pair_id": ["sr1"] * 4,
                "lr_pair_id": ["lr1"] * 4,
                "observed_mean_score": [1.0, 1.0, 1.0, 1.0],
                "null_mean": [0.5, 0.5, 0.5, 0.5],
                "pvalue": [0.01, 0.01, 0.5, 0.5],
            }
        )
        section_to_patient = {"S1": "P1", "S2": "P2", "S3": "P3", "S4": "P4"}
        result = aggregate_spatial_scores(scores, section_to_patient)
        row = result.iloc[0]
        assert row["n_sections_tested"] == 4
        assert row["n_sections_significant"] == 2
        assert row["spatial_specificity"] == 0.5

    def test_real_cross_patient_consistency_deduplicates_replicate_sections(self):
        # P1 has 2 replicate sections, only one significant -- should
        # still count as 1 real patient with a significant hit, not 2.
        scores = pd.DataFrame(
            {
                "section_id": ["S1a", "S1b", "S2"],
                "sender_receiver_pair_id": ["sr1"] * 3,
                "lr_pair_id": ["lr1"] * 3,
                "observed_mean_score": [1.0, 1.0, 1.0],
                "null_mean": [0.5, 0.5, 0.5],
                "pvalue": [0.01, 0.5, 0.5],
            }
        )
        section_to_patient = {"S1a": "P1", "S1b": "P1", "S2": "P2"}
        result = aggregate_spatial_scores(scores, section_to_patient)
        row = result.iloc[0]
        assert row["n_patients_tested"] == 2
        assert row["n_patients_with_significant_section"] == 1
        assert row["cross_patient_consistency"] == 0.5

    def test_real_effect_size_excludes_zero_null_or_observed_rows(self):
        scores = pd.DataFrame(
            {
                "section_id": ["S1", "S2", "S3"],
                "sender_receiver_pair_id": ["sr1"] * 3,
                "lr_pair_id": ["lr1"] * 3,
                "observed_mean_score": [2.0, 0.0, 4.0],
                "null_mean": [1.0, 0.0, 1.0],
                "pvalue": [0.01, np.nan, 0.01],
            }
        )
        section_to_patient = {"S1": "P1", "S2": "P2", "S3": "P3"}
        result = aggregate_spatial_scores(scores, section_to_patient)
        row = result.iloc[0]
        # log2(2/1)=1.0, log2(4/1)=2.0 -- real median of the two valid rows is 1.5.
        assert row["effect_size"] == 1.5
        assert row["n_sections_effect_size"] == 2


class TestComputeCompositePriority:
    def test_real_higher_criteria_values_rank_first(self):
        table = pd.DataFrame(
            {
                "pair_id": ["a", "b", "c"],
                "effect_size": [3.0, 1.0, 2.0],
                "spatial_specificity": [0.9, 0.1, 0.5],
            }
        )
        result = compute_composite_priority(table, ["effect_size", "spatial_specificity"])
        ordered = result.sort_values("priority_rank")["pair_id"].tolist()
        assert ordered == ["a", "c", "b"]

    def test_real_rank_1_is_assigned_to_top_row(self):
        table = pd.DataFrame({"pair_id": ["a", "b"], "effect_size": [5.0, 1.0]})
        result = compute_composite_priority(table, ["effect_size"])
        top = result[result["pair_id"] == "a"].iloc[0]
        assert top["priority_rank"] == 1
