"""Unit tests for xenium_tcr_ecology.qc.segmentation_quality (`04_quality_control/03_assess_segmentation_quality.py`)."""

from __future__ import annotations

from shapely.geometry import Polygon

from xenium_tcr_ecology.qc.segmentation_quality import check_cell_nucleus_pair


class TestCheckCellNucleusPair:
    def test_correct_containment_passes(self):
        cell = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        nucleus = Polygon([(3, 3), (7, 3), (7, 7), (3, 7)])
        result = check_cell_nucleus_pair(cell, nucleus)
        assert result["cell_valid"] and result["nucleus_valid"]
        assert result["nucleus_contained"]
        assert not result["nucleus_area_exceeds_cell"]

    def test_nucleus_outside_cell_is_flagged(self):
        cell = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        nucleus = Polygon([(20, 20), (24, 20), (24, 24), (20, 24)])  # entirely elsewhere
        result = check_cell_nucleus_pair(cell, nucleus)
        assert not result["nucleus_contained"]
        assert result["nucleus_overlap_fraction"] == 0.0

    def test_shared_boundary_edge_is_not_flagged(self):
        """Regression test: real Xenium cell boundaries are built by expanding
        outward from the nucleus boundary, so the two polygons legitimately
        share edges/vertices. A strict shapely `.contains()` check flags these
        as "not contained" purely from floating-point boundary-touching --
        confirmed against real data (~43% false-flag rate). A nucleus that
        shares an edge with its cell but is otherwise fully inside must pass."""
        cell = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        nucleus = Polygon([(0, 0), (10, 0), (10, 5), (0, 5)])  # shares the bottom edge with cell
        result = check_cell_nucleus_pair(cell, nucleus)
        assert result["nucleus_contained"]
        assert result["nucleus_overlap_fraction"] == 1.0

    def test_mostly_but_not_substantially_overlapping_nucleus_is_flagged(self):
        cell = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        # Nucleus mostly overlaps the cell but a large chunk sticks out -- a
        # genuine, non-trivial containment failure (overlap well below the
        # 0.95 threshold), unlike the shared-boundary sliver case above.
        nucleus = Polygon([(5, 5), (15, 5), (15, 15), (5, 15)])
        result = check_cell_nucleus_pair(cell, nucleus)
        assert not result["nucleus_contained"]
        assert result["nucleus_overlap_fraction"] < 0.95

    def test_nucleus_larger_than_cell_is_flagged(self):
        cell = Polygon([(0, 0), (5, 0), (5, 5), (0, 5)])
        nucleus = Polygon([(-5, -5), (15, -5), (15, 15), (-5, 15)])  # larger, engulfs cell
        result = check_cell_nucleus_pair(cell, nucleus)
        assert result["nucleus_area_exceeds_cell"]

    def test_invalid_polygon_detected(self):
        # A self-intersecting "bowtie" polygon is invalid.
        bowtie = Polygon([(0, 0), (10, 10), (10, 0), (0, 10)])
        cell = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        result = check_cell_nucleus_pair(cell, bowtie)
        assert not result["nucleus_valid"]
        # containment/area checks are skipped (default False) when either polygon is invalid
        assert result["nucleus_contained"] is False
