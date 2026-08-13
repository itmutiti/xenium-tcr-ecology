"""Unit tests for xenium_tcr_ecology.qc.transcript_spillover (`04_quality_control/04_estimate_transcript_spillover.py`)."""

from __future__ import annotations

import pandas as pd
import pytest
from shapely.geometry import Polygon

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.qc.transcript_spillover import compute_section_spillover_risk


def _square(x0: float, y0: float, side: float = 5.0) -> Polygon:
    return Polygon([(x0, y0), (x0 + side, y0), (x0 + side, y0 + side), (x0, y0 + side)])


class TestComputeSectionSpilloverRisk:
    def test_raises_on_mismatched_index(self):
        polys = pd.Series([_square(0, 0)], index=["a"])
        lineages = pd.Series(["T_cell"], index=["b"])
        with pytest.raises(PipelineError, match="same index"):
            compute_section_spillover_risk(polys, lineages)

    def test_empty_input_returns_empty_frame(self):
        result = compute_section_spillover_risk(pd.Series(dtype=object), pd.Series(dtype=object))
        assert len(result) == 0

    def test_isolated_cell_has_zero_risk(self):
        # Two far-apart same-type cells and one isolated cell -- the isolated
        # one has no neighbours within radius at all.
        polys = pd.Series(
            [_square(0, 0), _square(5, 0), _square(1000, 1000)], index=["a", "b", "isolated"]
        )
        lineages = pd.Series(["T_cell", "T_cell", "Myeloid"], index=["a", "b", "isolated"])
        result = compute_section_spillover_risk(polys, lineages, search_radius_um=10.0)
        assert result.loc["isolated", "n_neighbors_within_radius"] == 0
        assert result.loc["isolated", "spillover_risk_score"] == 0.0
        assert not result.loc["isolated", "is_boundary_adjacent_to_different_type"]

    def test_same_type_touching_neighbor_contributes_no_risk(self):
        polys = pd.Series([_square(0, 0), _square(5, 0)], index=["a", "b"])
        lineages = pd.Series(["T_cell", "T_cell"], index=["a", "b"])
        result = compute_section_spillover_risk(polys, lineages, search_radius_um=10.0)
        assert result.loc["a", "spillover_risk_score"] == 0.0
        assert result.loc["a", "n_different_type_neighbors_within_radius"] == 0
        assert not result.loc["a", "is_boundary_adjacent_to_different_type"]

    def test_touching_different_type_neighbor_flagged_and_scored(self):
        polys = pd.Series([_square(0, 0), _square(5, 0)], index=["a", "b"])
        lineages = pd.Series(["T_cell", "Myeloid"], index=["a", "b"])
        result = compute_section_spillover_risk(polys, lineages, search_radius_um=10.0)
        assert result.loc["a", "is_boundary_adjacent_to_different_type"]
        assert result.loc["a", "nearest_different_type_distance_um"] == 0.0
        # Sole neighbour, touching (distance 0) and different type: weight = 1 - 0/10 = 1.0.
        assert result.loc["a", "spillover_risk_score"] == pytest.approx(1.0)

    def test_risk_score_decays_with_distance(self):
        # Same different-type neighbour configuration, but pushed further
        # away (still within radius) -- risk score should be strictly lower.
        near = pd.Series([_square(0, 0), _square(5, 0)], index=["a", "b"])
        near_lineages = pd.Series(["T_cell", "Myeloid"], index=["a", "b"])
        near_result = compute_section_spillover_risk(near, near_lineages, search_radius_um=10.0)

        far = pd.Series([_square(0, 0), _square(8, 0)], index=["a", "b"])
        far_lineages = pd.Series(["T_cell", "Myeloid"], index=["a", "b"])
        far_result = compute_section_spillover_risk(far, far_lineages, search_radius_um=10.0)

        assert (
            far_result.loc["a", "spillover_risk_score"]
            < near_result.loc["a", "spillover_risk_score"]
        )
        assert not far_result.loc["a", "is_boundary_adjacent_to_different_type"]

    def test_neighbor_beyond_radius_is_ignored(self):
        polys = pd.Series([_square(0, 0), _square(50, 0)], index=["a", "b"])
        lineages = pd.Series(["T_cell", "Myeloid"], index=["a", "b"])
        result = compute_section_spillover_risk(polys, lineages, search_radius_um=10.0)
        assert result.loc["a", "n_neighbors_within_radius"] == 0
        assert result.loc["a", "spillover_risk_score"] == 0.0

    def test_multiple_different_type_neighbors_increase_risk(self):
        # 'a' is surrounded by two different-type touching neighbours instead
        # of one -- more different-type neighbours should not decrease risk.
        one_neighbor = pd.Series([_square(0, 0), _square(5, 0)], index=["a", "b"])
        one_lineages = pd.Series(["T_cell", "Myeloid"], index=["a", "b"])
        one_result = compute_section_spillover_risk(
            one_neighbor, one_lineages, search_radius_um=10.0
        )

        two_neighbors = pd.Series(
            [_square(0, 0), _square(5, 0), _square(-5, 0)], index=["a", "b", "c"]
        )
        two_lineages = pd.Series(["T_cell", "Myeloid", "B_cell"], index=["a", "b", "c"])
        two_result = compute_section_spillover_risk(
            two_neighbors, two_lineages, search_radius_um=10.0
        )

        assert (
            two_result.loc["a", "spillover_risk_score"]
            >= one_result.loc["a", "spillover_risk_score"]
        )
