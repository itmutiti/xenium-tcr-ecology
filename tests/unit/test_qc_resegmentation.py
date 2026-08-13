"""Unit tests for xenium_tcr_ecology.qc.resegmentation (`04_quality_control/05_resegment_reference_subset.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.qc.resegmentation import (
    MAX_EXPANSION_RADIUS_UM,
    reassign_transcripts_to_nearest_nucleus,
    summarize_reassignment_concordance,
)


class TestReassignTranscriptsToNearestNucleus:
    def test_assigns_to_nearest_nucleus(self):
        nucleus_ids = np.array(["a", "b"])
        nucleus_xy = np.array([[0.0, 0.0], [20.0, 0.0]])
        transcript_xy = np.array([[1.0, 0.0], [19.0, 0.0]])
        result = reassign_transcripts_to_nearest_nucleus(
            transcript_xy, nucleus_ids, nucleus_xy, max_radius_um=15.0
        )
        assert list(result) == ["a", "b"]

    def test_transcript_beyond_radius_is_background(self):
        nucleus_ids = np.array(["a"])
        nucleus_xy = np.array([[0.0, 0.0]])
        transcript_xy = np.array([[100.0, 0.0]])
        result = reassign_transcripts_to_nearest_nucleus(
            transcript_xy, nucleus_ids, nucleus_xy, max_radius_um=15.0
        )
        assert result[0] is None

    def test_transcript_within_radius_is_assigned(self):
        nucleus_ids = np.array(["a"])
        nucleus_xy = np.array([[0.0, 0.0]])
        transcript_xy = np.array([[10.0, 0.0]])
        result = reassign_transcripts_to_nearest_nucleus(
            transcript_xy, nucleus_ids, nucleus_xy, max_radius_um=15.0
        )
        assert result[0] == "a"

    def test_no_nuclei_returns_all_background(self):
        result = reassign_transcripts_to_nearest_nucleus(
            np.array([[0.0, 0.0]]), np.array([]), np.empty((0, 2)), max_radius_um=15.0
        )
        assert result[0] is None

    def test_no_transcripts_returns_empty(self):
        result = reassign_transcripts_to_nearest_nucleus(
            np.empty((0, 2)), np.array(["a"]), np.array([[0.0, 0.0]]), max_radius_um=15.0
        )
        assert len(result) == 0

    def test_default_radius_matches_documented_10x_default(self):
        assert MAX_EXPANSION_RADIUS_UM == 15.0


class TestSummarizeReassignmentConcordance:
    def test_fully_concordant(self):
        table = pd.DataFrame({"cell_id": ["x", "y"], "reassigned_nucleus_id": ["x", "y"]})
        result = summarize_reassignment_concordance(table)
        assert result["fraction_concordant_same_cell"] == 1.0
        assert result["fraction_primary_assigned_alt_background"] == 0.0

    def test_discordant_reassignment(self):
        table = pd.DataFrame({"cell_id": ["x", "y"], "reassigned_nucleus_id": ["z", "y"]})
        result = summarize_reassignment_concordance(table)
        assert result["fraction_concordant_same_cell"] == 0.5

    def test_primary_assigned_alternative_background(self):
        table = pd.DataFrame({"cell_id": ["x", "y"], "reassigned_nucleus_id": [None, "y"]})
        result = summarize_reassignment_concordance(table)
        assert result["fraction_primary_assigned_alt_background"] == 0.5

    def test_primary_unassigned_alternative_recovered(self):
        table = pd.DataFrame({"cell_id": ["UNASSIGNED", "y"], "reassigned_nucleus_id": ["x", "y"]})
        result = summarize_reassignment_concordance(table)
        assert result["fraction_primary_unassigned_alt_recovered"] == 0.5

    def test_both_unassigned(self):
        table = pd.DataFrame({"cell_id": ["UNASSIGNED"], "reassigned_nucleus_id": [None]})
        result = summarize_reassignment_concordance(table)
        assert result["fraction_both_unassigned"] == 1.0
        assert result["fraction_concordant_same_cell"] is None
