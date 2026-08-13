"""Unit tests for xenium_tcr_ecology.graphs.primary_graph (`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from xenium_tcr_ecology.graphs.primary_graph import build_weighted_graph


class TestBuildWeightedGraph:
    def test_edge_weight_is_real_euclidean_distance(self):
        coords = np.array([[0.0, 0.0], [3.0, 4.0]])
        adjacency = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        weighted = build_weighted_graph(adjacency, coords)
        assert weighted[0, 1] == 5.0  # 3-4-5 triangle
        assert weighted[1, 0] == 5.0

    def test_preserves_graph_topology(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [100.0, 100.0]])
        adjacency = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(3, 3))
        weighted = build_weighted_graph(adjacency, coords)
        assert weighted[0, 2] == 0  # no edge -> stays no edge
        assert weighted.nnz == adjacency.nnz

    def test_output_shape_matches_input(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0]])
        adjacency = sparse.csr_matrix(([1.0], ([0], [1])), shape=(2, 2))
        weighted = build_weighted_graph(adjacency, coords)
        assert weighted.shape == adjacency.shape
