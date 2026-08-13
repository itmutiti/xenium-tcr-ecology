"""Unit tests for xenium_tcr_ecology.niches.local_composition (`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`)."""

from __future__ import annotations

import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.niches.local_composition import compute_composition_vectors


class TestComputeCompositionVectors:
    def test_pure_neighbourhood_gives_100_percent_that_lineage(self):
        # Cell 0 (T_cell) surrounded by 3 B_cell neighbours.
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0, 1.0, 1.0], ([0, 0, 0, 1, 2, 3], [1, 2, 3, 0, 0, 0])), shape=(4, 4)
        )
        lineage = pd.Series(
            ["T_cell", "B_cell", "B_cell", "B_cell"], index=["c0", "c1", "c2", "c3"]
        )
        result = compute_composition_vectors(graph, lineage, ["B_cell", "T_cell"])
        assert result.loc["c0", "B_cell"] == 1.0
        assert result.loc["c0", "T_cell"] == 0.0

    def test_own_lineage_is_excluded_from_own_composition(self):
        # Cell 0 (B_cell) has one B_cell neighbour and one T_cell neighbour
        # -- composition should be 50/50, not counting cell 0 itself.
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 0, 1, 2], [1, 2, 0, 0])), shape=(3, 3)
        )
        lineage = pd.Series(["B_cell", "B_cell", "T_cell"], index=["c0", "c1", "c2"])
        result = compute_composition_vectors(graph, lineage, ["B_cell", "T_cell"])
        assert result.loc["c0", "B_cell"] == 0.5
        assert result.loc["c0", "T_cell"] == 0.5

    def test_zero_degree_cell_gets_nan_not_zero(self):
        graph = sparse.csr_matrix((2, 2))
        lineage = pd.Series(["B_cell", "T_cell"], index=["c0", "c1"])
        result = compute_composition_vectors(graph, lineage, ["B_cell", "T_cell"])
        assert result.loc["c0"].isna().all()

    def test_rows_sum_to_one_for_connected_cells(self):
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 0, 1, 2], [1, 2, 0, 0])), shape=(3, 3)
        )
        lineage = pd.Series(["B_cell", "T_cell", "B_cell"], index=["c0", "c1", "c2"])
        result = compute_composition_vectors(graph, lineage, ["B_cell", "T_cell"])
        assert result.loc["c0"].sum() == 1.0
