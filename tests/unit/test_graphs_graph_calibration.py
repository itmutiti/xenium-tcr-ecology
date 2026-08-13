"""Unit tests for xenium_tcr_ecology.graphs.graph_calibration (`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`)."""

from __future__ import annotations

from scipy import sparse

from xenium_tcr_ecology.graphs.graph_calibration import (
    compute_largest_component_fraction,
    select_calibrated_parameter,
)


class TestComputeLargestComponentFraction:
    def test_fully_connected_graph_returns_one(self):
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 1, 1, 2], [1, 0, 2, 1])), shape=(3, 3)
        )
        assert compute_largest_component_fraction(graph) == 1.0

    def test_split_graph_returns_the_larger_fraction(self):
        # Two components: {0,1} and {2,3} -- each half the graph.
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 2], [1, 3])), shape=(4, 4))
        assert compute_largest_component_fraction(graph) == 0.5

    def test_empty_graph_returns_one(self):
        graph = sparse.csr_matrix((0, 0))
        assert compute_largest_component_fraction(graph) == 1.0


class TestSelectCalibratedParameter:
    def test_picks_smallest_scale_among_those_clearing_threshold(self):
        connectivity = {"a": 0.99, "b": 0.97, "c": 0.5}
        scale = {"a": 30.0, "b": 15.0, "c": 50.0}
        result = select_calibrated_parameter(connectivity, scale, threshold=0.95)
        assert result == "b"  # smaller scale, still clears 0.95

    def test_falls_back_to_most_connected_when_none_clear_threshold(self):
        connectivity = {"a": 0.5, "b": 0.7, "c": 0.3}
        scale = {"a": 15.0, "b": 30.0, "c": 50.0}
        result = select_calibrated_parameter(connectivity, scale, threshold=0.95)
        assert result == "b"  # highest connectivity, even though it doesn't clear the bar

    def test_single_candidate_clearing_threshold_is_selected(self):
        connectivity = {"a": 0.5, "b": 0.99}
        scale = {"a": 15.0, "b": 30.0}
        result = select_calibrated_parameter(connectivity, scale, threshold=0.95)
        assert result == "b"
