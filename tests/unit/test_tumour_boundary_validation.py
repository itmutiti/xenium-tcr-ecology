"""Unit tests for xenium_tcr_ecology.tumour.boundary_validation (`07_tumour_epithelium_characterisation/06_validate_boundaries_against_morphology.py`)."""

from __future__ import annotations

from shapely.geometry import Polygon

from xenium_tcr_ecology.tumour.boundary_validation import (
    REVIEW_LOG_COLUMNS,
    anonymize_panel_id,
    sample_boundary_points,
)


class TestAnonymizePanelId:
    def test_deterministic(self):
        assert anonymize_panel_id("P01_run1", 0) == anonymize_panel_id("P01_run1", 0)

    def test_different_index_gives_different_id(self):
        assert anonymize_panel_id("P01_run1", 0) != anonymize_panel_id("P01_run1", 1)

    def test_different_section_gives_different_id(self):
        assert anonymize_panel_id("P01_run1", 0) != anonymize_panel_id("P09_run1", 0)


class TestSampleBoundaryPoints:
    def test_returns_requested_number_of_points(self):
        square = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        points = sample_boundary_points(square, n_points=4)
        assert len(points) == 4

    def test_points_lie_on_the_boundary(self):
        from shapely.geometry import Point

        square = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        points = sample_boundary_points(square, n_points=4)
        boundary = square.boundary
        for x, y in points:
            assert boundary.distance(Point(x, y)) < 1e-6

    def test_zero_length_boundary_returns_empty(self):
        from shapely.geometry import Point

        degenerate = Point(0, 0).buffer(0)  # empty polygon
        points = sample_boundary_points(degenerate, n_points=4)
        assert points == []


def test_review_log_columns_include_panel_id_and_no_prefilled_agreement():
    assert "panel_id" in REVIEW_LOG_COLUMNS
    assert "reviewer_agrees_with_boundary" in REVIEW_LOG_COLUMNS
