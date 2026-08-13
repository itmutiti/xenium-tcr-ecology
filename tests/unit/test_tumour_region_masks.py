"""Unit tests for xenium_tcr_ecology.tumour.region_masks (`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.tumour.region_masks import (
    compute_smoothed_calls,
    filter_small_regions,
    label_connected_regions,
)


class TestComputeSmoothedCalls:
    def test_isolated_positive_call_surrounded_by_negatives_is_smoothed_away(self):
        # A tight cluster of 6 cells, 5 negative + 1 positive at the centre
        # -- the positive cell's neighbours are all negative, so majority
        # vote should smooth it to False.
        x = np.array([0.0, 1.0, -1.0, 0.0, 0.0, 0.5])
        y = np.array([0.0, 0.0, 0.0, 1.0, -1.0, 0.5])
        raw_calls = np.array([True, False, False, False, False, False])
        result = compute_smoothed_calls(x, y, raw_calls, k=4)
        assert not result[0]

    def test_majority_positive_neighbourhood_stays_positive(self):
        x = np.array([0.0, 1.0, -1.0, 0.0, 0.0])
        y = np.array([0.0, 0.0, 0.0, 1.0, -1.0])
        raw_calls = np.array([True, True, True, True, False])
        result = compute_smoothed_calls(x, y, raw_calls, k=4)
        assert result[0]

    def test_empty_input_returns_empty(self):
        result = compute_smoothed_calls(np.array([]), np.array([]), np.array([], dtype=bool))
        assert len(result) == 0


class TestLabelConnectedRegions:
    def test_two_separate_clusters_get_different_region_ids(self):
        x = np.array([0.0, 1.0, 100.0, 101.0])
        y = np.array([0.0, 0.0, 0.0, 0.0])
        is_malignant = np.array([True, True, True, True])
        result = label_connected_regions(x, y, is_malignant, radius_um=5.0)
        assert result[0] == result[1]
        assert result[2] == result[3]
        assert result[0] != result[2]

    def test_non_malignant_cells_get_no_region(self):
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 0.0])
        is_malignant = np.array([True, False])
        result = label_connected_regions(x, y, is_malignant, radius_um=5.0)
        assert result[1] == -1
        assert result[0] != -1

    def test_no_malignant_cells_returns_all_unassigned(self):
        x = np.array([0.0, 1.0])
        y = np.array([0.0, 0.0])
        is_malignant = np.array([False, False])
        result = label_connected_regions(x, y, is_malignant)
        assert (result == -1).all()


class TestFilterSmallRegions:
    def test_removes_regions_below_min_size(self):
        region_id = np.array([0, 0, 0, 0, 0, 1, -1])  # region 0 has 5 cells, region 1 has 1
        result = filter_small_regions(region_id, min_size=5)
        assert (result[:5] == 0).all()
        assert result[5] == -1
        assert result[6] == -1

    def test_no_valid_regions_is_a_noop(self):
        region_id = np.array([-1, -1, -1])
        result = filter_small_regions(region_id, min_size=5)
        assert (result == -1).all()

    def test_region_at_exactly_min_size_is_kept(self):
        region_id = np.array([0, 0, 0])
        result = filter_small_regions(region_id, min_size=3)
        assert (result == 0).all()
