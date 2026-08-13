"""Unit tests for xenium_tcr_ecology.graphs.graph_pruning (`09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from xenium_tcr_ecology.graphs.graph_pruning import (
    compute_max_edge_length,
    modified_z_scores,
    prune_long_edges,
)


class TestModifiedZScores:
    def test_flags_real_outlier(self):
        values = np.array([10.0, 10.5, 9.8, 10.2, 9.9, 50.0])
        scores = modified_z_scores(values)
        assert abs(scores[-1]) > 3.5
        assert all(abs(s) < 3.5 for s in scores[:-1])

    def test_zero_mad_returns_zeros(self):
        values = np.array([5.0, 5.0, 5.0, 5.0])
        scores = modified_z_scores(values)
        assert (scores == 0).all()


class TestComputeMaxEdgeLength:
    def test_empty_input_returns_infinity(self):
        assert compute_max_edge_length(np.array([])) == np.inf

    def test_threshold_is_grounded_in_the_real_distribution(self):
        # A continuously-varying cluster around 15 (real MAD > 0, unlike a
        # handful of repeated discrete values) with one dramatic outlier
        # -- the threshold should sit well above the cluster but not be
        # driven up close to the outlier itself (robust to the outlier,
        # the whole point of using MAD rather than mean+std).
        rng = np.random.default_rng(0)
        cluster = rng.normal(15.0, 3.0, size=200)
        lengths = np.concatenate([cluster, [1000.0]])
        threshold = compute_max_edge_length(lengths)
        assert 15 < threshold < 100

    def test_zero_mad_falls_back_to_max(self):
        lengths = np.array([10.0, 10.0, 10.0, 500.0])
        # median=10, MAD=0 (three identical values dominate) -- degenerate
        # case must not divide by zero or return something nonsensical.
        threshold = compute_max_edge_length(lengths)
        assert threshold == 500.0


class TestPruneLongEdges:
    def test_removes_edges_above_threshold(self):
        coords = np.array([[0.0, 0.0], [10.0, 0.0], [1000.0, 0.0]])
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 1, 1, 2], [1, 0, 2, 1])), shape=(3, 3)
        )
        pruned = prune_long_edges(graph, coords, max_edge_length_um=50.0)
        assert pruned[0, 1] != 0  # short edge kept
        assert pruned[1, 2] == 0  # long edge (990um) removed

    def test_keeps_all_edges_when_threshold_is_generous(self):
        coords = np.array([[0.0, 0.0], [10.0, 0.0]])
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        pruned = prune_long_edges(graph, coords, max_edge_length_um=1000.0)
        assert pruned.nnz == graph.nnz

    def test_output_shape_matches_input(self):
        coords = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0]])
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(3, 3))
        pruned = prune_long_edges(graph, coords, max_edge_length_um=5.0)
        assert pruned.shape == graph.shape
