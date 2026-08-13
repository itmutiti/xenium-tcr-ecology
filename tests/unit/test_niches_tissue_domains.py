"""Unit tests for xenium_tcr_ecology.niches.tissue_domains (`10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from xenium_tcr_ecology.niches.tissue_domains import (
    find_contiguous_domains,
    smooth_labels_by_majority_vote,
)


def _path_graph(n: int) -> sparse.csr_matrix:
    rows = list(range(n - 1)) + list(range(1, n))
    cols = list(range(1, n)) + list(range(n - 1))
    return sparse.csr_matrix(([1.0] * len(rows), (rows, cols)), shape=(n, n))


class TestSmoothLabelsByMajorityVote:
    def test_isolated_salt_cell_is_relabelled_to_neighbour_majority(self):
        # 5-cell path graph, all label 0 except a single "salt" cell (index 2, label 1)
        # surrounded on both sides by label-0 neighbours.
        graph = _path_graph(5)
        labels = np.array([0, 0, 1, 0, 0])
        smoothed = smooth_labels_by_majority_vote(graph, labels)
        assert smoothed[2] == 0

    def test_consistent_block_is_unchanged(self):
        graph = _path_graph(5)
        labels = np.array([0, 0, 0, 0, 0])
        smoothed = smooth_labels_by_majority_vote(graph, labels)
        assert np.array_equal(smoothed, labels)

    def test_tie_among_neighbours_is_broken_by_the_self_vote(self):
        # Star graph: centre cell (0, label 1) connects to two label-0 and
        # two label-1 neighbours -- neighbours alone split 2 vs 2, but the
        # centre's own self-vote (included in the vote count) tips the
        # count for its own label (3) above the alternative (2).
        rows = [0, 0, 0, 0]
        cols = [1, 2, 3, 4]
        graph = sparse.csr_matrix(([1.0] * 8, (rows + cols, cols + rows)), shape=(5, 5))
        labels = np.array(
            [1, 0, 0, 1, 1]
        )  # centre=1; neighbours: 0,0,1,1 -> 2 vs 2 among neighbours
        smoothed = smooth_labels_by_majority_vote(graph, labels)
        assert smoothed[0] == 1  # own label, wins via self-vote, kept

    def test_zero_degree_cell_keeps_own_label(self):
        graph = sparse.csr_matrix((2, 2))
        labels = np.array([0, 1])
        smoothed = smooth_labels_by_majority_vote(graph, labels)
        assert np.array_equal(smoothed, labels)


class TestFindContiguousDomains:
    def test_unconnected_cells_with_the_same_label_are_not_merged(self):
        # Only cells 0-1 are graph-connected (both label 0); cells 2 and 3
        # have no edges to anything, including each other. Sharing a label
        # must not merge cells that have no real spatial adjacency.
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(4, 4))
        labels = np.array([0, 0, 0, 0])
        domains = find_contiguous_domains(graph, labels)
        assert domains[0] == domains[1]  # real graph edge, same label -> same domain
        assert (
            domains[2] != domains[3]
        )  # no edge at all -> distinct singleton domains despite same label

    def test_different_label_neighbours_are_split_into_separate_domains(self):
        graph = _path_graph(4)
        labels = np.array([0, 0, 1, 1])
        domains = find_contiguous_domains(graph, labels)
        assert domains[0] == domains[1]
        assert domains[2] == domains[3]
        assert domains[0] != domains[2]

    def test_rare_single_cell_niche_keeps_its_own_domain_not_merged(self):
        # 3-cell path, middle cell has a different label from both neighbours
        # (which happen to share a label) -- middle cell must be its own
        # domain, not merged into either neighbour's domain.
        graph = _path_graph(3)
        labels = np.array([0, 1, 0])
        domains = find_contiguous_domains(graph, labels)
        assert domains[1] != domains[0]
        assert domains[1] != domains[2]
        assert domains[0] != domains[2]  # the two label-0 cells aren't directly connected either
