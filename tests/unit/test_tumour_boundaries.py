"""Unit tests for xenium_tcr_ecology.tumour.tumour_boundaries (`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`)."""

from __future__ import annotations

import numpy as np
import pytest

from xenium_tcr_ecology.tumour.tumour_boundaries import (
    build_tumour_mask_polygon,
    compute_signed_distances,
)


class TestBuildTumourMaskPolygon:
    def test_empty_input_returns_none(self):
        result = build_tumour_mask_polygon(np.array([]), np.array([]))
        assert result is None

    def test_single_point_produces_a_disk(self):
        result = build_tumour_mask_polygon(np.array([0.0]), np.array([0.0]), buffer_um=5.0)
        assert result is not None
        assert result.area == pytest.approx(np.pi * 5.0**2, rel=0.01)

    def test_nearby_points_merge_into_one_polygon(self):
        result = build_tumour_mask_polygon(
            np.array([0.0, 1.0]), np.array([0.0, 0.0]), buffer_um=5.0
        )
        # Two overlapping disks merge -- area is less than the sum of two
        # separate disks (no double-counted overlap), but the union is a
        # single connected polygon, not two.
        assert result.geom_type == "Polygon"
        assert result.area < 2 * np.pi * 5.0**2


class TestComputeSignedDistances:
    def test_no_polygon_returns_null_and_not_inside(self):
        result = compute_signed_distances(np.array([0.0, 1.0]), np.array([0.0, 1.0]), None)
        assert result["distance_to_tumour_boundary_um"].isna().all()
        assert not result["is_inside_tumour_region"].any()

    def test_point_inside_gets_negative_signed_distance(self):
        polygon = build_tumour_mask_polygon(np.array([0.0]), np.array([0.0]), buffer_um=10.0)
        result = compute_signed_distances(np.array([0.0]), np.array([0.0]), polygon)
        assert result.loc[0, "is_inside_tumour_region"]
        assert result.loc[0, "signed_distance_to_tumour_boundary_um"] < 0
        assert result.loc[0, "distance_to_tumour_boundary_um"] == pytest.approx(10.0, rel=0.01)

    def test_point_outside_gets_positive_signed_distance(self):
        polygon = build_tumour_mask_polygon(np.array([0.0]), np.array([0.0]), buffer_um=5.0)
        result = compute_signed_distances(np.array([100.0]), np.array([0.0]), polygon)
        assert not result.loc[0, "is_inside_tumour_region"]
        assert result.loc[0, "signed_distance_to_tumour_boundary_um"] > 0
        assert result.loc[0, "distance_to_tumour_boundary_um"] == pytest.approx(95.0, rel=0.01)

    def test_unsigned_distance_is_always_the_absolute_value_of_signed(self):
        polygon = build_tumour_mask_polygon(np.array([0.0]), np.array([0.0]), buffer_um=8.0)
        result = compute_signed_distances(np.array([0.0, 50.0]), np.array([0.0, 0.0]), polygon)
        assert (
            result["distance_to_tumour_boundary_um"]
            == result["signed_distance_to_tumour_boundary_um"].abs()
        ).all()
