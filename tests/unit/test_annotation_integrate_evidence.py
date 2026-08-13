"""Unit tests for xenium_tcr_ecology.annotation.integrate_evidence (`06_cell_type_annotation/06_integrate_annotation_evidence.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.integrate_evidence import (
    compute_cluster_consensus,
    compute_marker_margin,
    compute_spatial_consistency,
    integrate_evidence,
)


class TestComputeMarkerMargin:
    def test_larger_gap_gives_higher_margin(self):
        scores = pd.DataFrame(
            {
                "A_lineage_score": [1.0, 1.0],
                "B_lineage_score": [0.9, 0.1],  # cell 0: small gap; cell 1: large gap
            },
            index=["c0", "c1"],
        )
        result = compute_marker_margin(scores)
        assert result["c1"] > result["c0"]

    def test_output_is_in_zero_one_range(self):
        rng = np.random.default_rng(0)
        scores = pd.DataFrame(
            rng.normal(size=(50, 4)), columns=[f"L{i}_lineage_score" for i in range(4)]
        )
        result = compute_marker_margin(scores)
        assert result.between(0, 1).all()

    def test_constant_margin_returns_a_finite_tied_value_not_nan(self):
        # Percentile rank (not min-max) is used precisely because it is
        # robust to skew in the raw margin distribution -- see the
        # function's own docstring/comment for the real-data motivation.
        # Tied inputs correctly get the same (non-zero) percentile rank,
        # not 0: a value of 0 would incorrectly signal "the lowest-
        # confidence cells in the population" when they are actually all
        # equally-ranked, not distinguishably low.
        scores = pd.DataFrame({"A_lineage_score": [1.0, 1.0], "B_lineage_score": [0.5, 0.5]})
        result = compute_marker_margin(scores)
        assert not result.isna().any()
        assert result.iloc[0] == result.iloc[1]
        assert result.between(0, 1).all()

    def test_uses_percentile_rank_not_min_max(self):
        # A single extreme outlier must not compress everyone else's
        # margin toward 0 (the min-max failure mode found on real data).
        scores = pd.DataFrame(
            {
                "A_lineage_score": [0.5, 0.5, 0.5, 100.0],
                "B_lineage_score": [0.4, 0.3, 0.2, 0.0],  # margins: 0.1, 0.2, 0.3, 100.0
            }
        )
        result = compute_marker_margin(scores)
        # Under min-max, the first three (raw margins 0.1/0.2/0.3) would
        # all be compressed below 0.003 by the outlier's raw margin of 100.
        # Percentile rank must not do this -- they should occupy a
        # reasonable, spread-out portion of [0, 1].
        assert result.iloc[0] > 0.1


class TestComputeClusterConsensus:
    def test_uniform_cluster_gives_consensus_of_one(self):
        clustering = pd.DataFrame(
            {"joint_leiden_res0.5": ["0", "0", "0"]}, index=["c0", "c1", "c2"]
        )
        lineage = pd.Series(["T_cell", "T_cell", "T_cell"], index=["c0", "c1", "c2"])
        result = compute_cluster_consensus(clustering, lineage)
        assert (result == 1.0).all()

    def test_minority_cell_gets_lower_consensus(self):
        clustering = pd.DataFrame(
            {"joint_leiden_res0.5": ["0", "0", "0", "0"]}, index=["c0", "c1", "c2", "c3"]
        )
        lineage = pd.Series(
            ["T_cell", "T_cell", "T_cell", "B_cell"], index=["c0", "c1", "c2", "c3"]
        )
        result = compute_cluster_consensus(clustering, lineage)
        assert result["c3"] < result["c0"]
        assert result["c0"] == 0.75

    def test_raises_if_resolution_column_missing(self):
        clustering = pd.DataFrame({"other_col": ["0"]}, index=["c0"])
        lineage = pd.Series(["T_cell"], index=["c0"])
        with pytest.raises(PipelineError, match="joint_leiden_res0.5"):
            compute_cluster_consensus(clustering, lineage)


class TestComputeSpatialConsistency:
    def test_isolated_cell_in_a_homogeneous_neighbourhood_scores_low(self):
        # A tight cluster of "T_cell" cells with one "B_cell" outlier
        # embedded right in the middle of it.
        rng = np.random.default_rng(0)
        n = 30
        x = rng.normal(0, 1, size=n)
        y = rng.normal(0, 1, size=n)
        labels = ["T_cell"] * n
        labels[0] = "B_cell"  # the cell at the cluster's own coordinates
        x[0], y[0] = 0.0, 0.0  # place it at the cluster centre, surrounded by T_cell neighbours

        section_ids = pd.Series(["S1"] * n, index=[f"c{i}" for i in range(n)])
        lineage = pd.Series(labels, index=section_ids.index)
        result = compute_spatial_consistency(section_ids, x, y, lineage, k=5)
        assert result["c0"] < 0.5

    def test_two_spatially_separated_homogeneous_clusters_score_high(self):
        n_per_cluster = 20
        x = np.concatenate(
            [
                np.random.default_rng(1).normal(0, 0.5, n_per_cluster),
                np.random.default_rng(2).normal(100, 0.5, n_per_cluster),
            ]
        )
        y = np.zeros(n_per_cluster * 2)
        labels = ["A"] * n_per_cluster + ["B"] * n_per_cluster
        section_ids = pd.Series(["S1"] * len(x), index=[f"c{i}" for i in range(len(x))])
        lineage = pd.Series(labels, index=section_ids.index)
        result = compute_spatial_consistency(section_ids, x, y, lineage, k=5)
        assert (result > 0.9).all()

    def test_sections_are_handled_independently(self):
        # Two sections, each internally consistent but with opposite
        # labels -- must not leak neighbours across sections.
        section_ids = pd.Series(["S1", "S1", "S2", "S2"], index=["c0", "c1", "c2", "c3"])
        x = np.array([0.0, 0.1, 0.0, 0.1])
        y = np.array([0.0, 0.1, 0.0, 0.1])
        lineage = pd.Series(["A", "A", "B", "B"], index=section_ids.index)
        result = compute_spatial_consistency(section_ids, x, y, lineage, k=1)
        assert (result == 1.0).all()


class TestIntegrateEvidence:
    def test_produces_expected_columns_and_no_nan_confidence(self):
        n = 20
        lineage_scores = pd.DataFrame(
            {"A_lineage_score": np.linspace(0, 1, n), "B_lineage_score": np.linspace(1, 0, n)},
            index=[f"c{i}" for i in range(n)],
        )
        clustering = pd.DataFrame({"joint_leiden_res0.5": ["0"] * n}, index=lineage_scores.index)
        section_ids = pd.Series(["S1"] * n, index=lineage_scores.index)
        x = np.linspace(0, 10, n)
        y = np.zeros(n)

        result = integrate_evidence(lineage_scores, clustering, section_ids, x, y)

        for col in [
            "final_lineage",
            "marker_margin",
            "cluster_consensus",
            "spatial_consistency",
            "confidence",
            "is_ambiguous",
            "final_substate",
        ]:
            assert col in result.columns
        assert not result["confidence"].isna().any()

    def test_ambiguous_cells_have_no_substate_even_if_one_was_provided(self):
        n = 4
        lineage_scores = pd.DataFrame(
            {
                "A_lineage_score": [0.0] * n,
                "B_lineage_score": [0.0] * n,
            },  # zero margin -> low confidence
            index=[f"c{i}" for i in range(n)],
        )
        clustering = pd.DataFrame(
            {"joint_leiden_res0.5": ["0", "1", "0", "1"]}, index=lineage_scores.index
        )
        section_ids = pd.Series(["S1"] * n, index=lineage_scores.index)
        x = np.array([0.0, 100.0, 0.1, 100.1])
        y = np.array([0.0, 0.0, 0.1, 0.1])
        substates = pd.DataFrame({"t_cell_state": ["Treg"] * n}, index=lineage_scores.index)

        result = integrate_evidence(
            lineage_scores, clustering, section_ids, x, y, t_cell_states=substates
        )
        assert result.loc[result["is_ambiguous"], "final_substate"].isna().all()
