"""Unit tests for xenium_tcr_ecology.graphs.tumour_tcell_bipartite (`09_spatial_graph_construction_and_calibration/04_construct_tumour_tcell_bipartite_graph.py`)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from xenium_tcr_ecology.graphs.tumour_tcell_bipartite import extract_bipartite_subgraph


class TestExtractBipartiteSubgraph:
    def test_keeps_cross_type_edge(self):
        # 0=malignant, 1=T cell
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        is_malignant = np.array([True, False])
        is_tcell = np.array([False, True])
        result = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
        assert result[0, 1] != 0
        assert result[1, 0] != 0

    def test_drops_malignant_malignant_edge(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        is_malignant = np.array([True, True])
        is_tcell = np.array([False, False])
        result = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
        assert result.nnz == 0

    def test_drops_tcell_tcell_edge(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        is_malignant = np.array([False, False])
        is_tcell = np.array([True, True])
        result = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
        assert result.nnz == 0

    def test_drops_edge_to_a_third_cell_type(self):
        # 0=malignant, 1=T cell, 2=fibroblast (neither malignant nor T cell)
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 1, 0, 2], [1, 0, 2, 0])), shape=(3, 3)
        )
        is_malignant = np.array([True, False, False])
        is_tcell = np.array([False, True, False])
        result = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
        assert result[0, 1] != 0  # malignant-Tcell kept
        assert result[0, 2] == 0  # malignant-fibroblast dropped

    def test_preserves_edge_weight(self):
        graph = sparse.csr_matrix(([7.5, 7.5], ([0, 1], [1, 0])), shape=(2, 2))
        is_malignant = np.array([True, False])
        is_tcell = np.array([False, True])
        result = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
        assert result[0, 1] == 7.5
