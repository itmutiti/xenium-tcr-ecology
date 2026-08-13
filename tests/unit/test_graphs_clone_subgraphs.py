"""Unit tests for xenium_tcr_ecology.graphs.clone_subgraphs (`09_spatial_graph_construction_and_calibration/05_construct_clone_induced_subgraphs.py`)."""

from __future__ import annotations

from scipy import sparse

from xenium_tcr_ecology.graphs.clone_subgraphs import compute_clone_shells


def _chain_graph(n):
    """A simple path graph 0-1-2-...-(n-1)."""
    rows = list(range(n - 1)) + list(range(1, n))
    cols = list(range(1, n)) + list(range(n - 1))
    data = [1.0] * len(rows)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


class TestComputeCloneShells:
    def test_clone_members_are_shell_zero(self):
        graph = _chain_graph(5)
        shells = compute_clone_shells({0}, graph, n_shells=3)
        assert shells[0] == 0

    def test_direct_neighbor_is_shell_one(self):
        graph = _chain_graph(5)
        shells = compute_clone_shells({0}, graph, n_shells=3)
        assert shells[1] == 1

    def test_two_hops_away_is_shell_two(self):
        graph = _chain_graph(5)
        shells = compute_clone_shells({0}, graph, n_shells=3)
        assert shells[2] == 2

    def test_cells_beyond_n_shells_are_absent(self):
        graph = _chain_graph(10)
        shells = compute_clone_shells({0}, graph, n_shells=2)
        assert 3 not in shells
        assert 2 in shells

    def test_multiple_clone_members_take_shortest_path(self):
        # Clone members at both ends of a 5-node chain -- the middle node
        # is 2 hops from either end, not further.
        graph = _chain_graph(5)
        shells = compute_clone_shells({0, 4}, graph, n_shells=3)
        assert shells[2] == 2
        assert shells[0] == 0
        assert shells[4] == 0

    def test_disconnected_cell_is_absent(self):
        # Node 5 has no edges at all.
        graph = _chain_graph(5)
        graph.resize((6, 6))
        shells = compute_clone_shells({0}, graph, n_shells=3)
        assert 5 not in shells
