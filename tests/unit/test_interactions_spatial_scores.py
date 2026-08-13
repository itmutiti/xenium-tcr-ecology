"""Unit tests for xenium_tcr_ecology.interactions.spatial_scores (`14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.interactions.spatial_scores import (
    compute_degree_preserving_lr_pvalue,
    compute_edge_interaction_scores,
    find_program_overlap_combinations,
)


class TestComputeEdgeInteractionScores:
    def test_real_edge_score_is_ligand_times_receptor(self):
        # Cell 0 (sender) -- Cell 1 (receiver), symmetric storage.
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        sender_mask = np.array([True, False])
        receiver_mask = np.array([False, True])
        ligand_expr = np.array([2.0, 0.0])
        receptor_expr = np.array([0.0, 3.0])
        result = compute_edge_interaction_scores(
            graph, sender_mask, receiver_mask, ligand_expr, receptor_expr
        )
        assert list(result) == [6.0]

    def test_no_sender_receiver_edges_gives_empty_array(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        sender_mask = np.array([True, True])
        receiver_mask = np.array([False, False])
        result = compute_edge_interaction_scores(
            graph, sender_mask, receiver_mask, np.array([1.0, 1.0]), np.array([1.0, 1.0])
        )
        assert len(result) == 0

    def test_both_edge_orientations_counted_from_symmetric_storage(self):
        # 3 cells: 0 (sender), 1 (receiver), 2 (sender) -- edges (0,1) and (2,1).
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 1, 2, 1], [1, 0, 1, 2])), shape=(3, 3)
        )
        sender_mask = np.array([True, False, True])
        receiver_mask = np.array([False, True, False])
        ligand_expr = np.array([2.0, 0.0, 3.0])
        receptor_expr = np.array([0.0, 5.0, 0.0])
        result = compute_edge_interaction_scores(
            graph, sender_mask, receiver_mask, ligand_expr, receptor_expr
        )
        assert sorted(result) == [10.0, 15.0]


class TestComputeDegreePreservingLrPvalue:
    def test_real_elevated_signal_gives_small_pvalue(self):
        rng_setup = np.random.default_rng(0)
        n = 200
        # A ring-like graph so every cell has degree 2 (uniform degree strata).
        rows = list(range(n))
        cols = [(i + 1) % n for i in range(n)]
        graph = sparse.csr_matrix(([1.0] * len(rows), (rows, cols)), shape=(n, n))
        graph = ((graph + graph.T) > 0).astype(float).tocsr()

        sender_mask = np.zeros(n, dtype=bool)
        receiver_mask = np.zeros(n, dtype=bool)
        sender_mask[:20] = True
        receiver_mask[1:21] = True  # adjacent to senders in the ring -> real elevated signal
        ligand_expr = np.where(sender_mask, 5.0, 0.1)
        receptor_expr = np.where(receiver_mask, 5.0, 0.1)
        degree_strata = np.zeros(n, dtype=int)  # uniform degree -> single stratum

        result = compute_degree_preserving_lr_pvalue(
            graph,
            sender_mask,
            receiver_mask,
            ligand_expr,
            receptor_expr,
            degree_strata,
            rng_setup,
            n_permutations=100,
        )
        assert result["n_edges"] > 0
        assert result["observed_mean_score"] > result["null_mean"]

    def test_no_real_edges_gives_nan_pvalue(self):
        graph = sparse.csr_matrix((4, 4))
        sender_mask = np.array([True, False, False, False])
        receiver_mask = np.array([False, True, False, False])
        rng_setup = np.random.default_rng(1)
        degree_strata = np.zeros(4, dtype=int)
        result = compute_degree_preserving_lr_pvalue(
            graph,
            sender_mask,
            receiver_mask,
            np.ones(4),
            np.ones(4),
            degree_strata,
            rng_setup,
            n_permutations=10,
        )
        assert np.isnan(result["pvalue"])


class TestFindProgramOverlapCombinations:
    def test_real_overlap_is_found(self):
        sender_receiver_pairs = [
            {"pair_id": "sr1", "sender": "A", "receiver": "B", "relevant_programs": ["exhaustion"]}
        ]
        lr_pairs = pd.DataFrame(
            [
                {
                    "pair_id": "lr1",
                    "ligand": "L",
                    "receptor": "R",
                    "programs": "exhaustion",
                    "pair_complete": True,
                }
            ]
        )
        result = find_program_overlap_combinations(sender_receiver_pairs, lr_pairs)
        assert len(result) == 1
        assert result[0]["lr_pair_id"] == "lr1"

    def test_no_overlap_is_excluded(self):
        sender_receiver_pairs = [
            {"pair_id": "sr1", "sender": "A", "receiver": "B", "relevant_programs": ["activation"]}
        ]
        lr_pairs = pd.DataFrame(
            [
                {
                    "pair_id": "lr1",
                    "ligand": "L",
                    "receptor": "R",
                    "programs": "exhaustion",
                    "pair_complete": True,
                }
            ]
        )
        result = find_program_overlap_combinations(sender_receiver_pairs, lr_pairs)
        assert result == []

    def test_incomplete_pairs_are_excluded_even_with_program_overlap(self):
        sender_receiver_pairs = [
            {"pair_id": "sr1", "sender": "A", "receiver": "B", "relevant_programs": ["exhaustion"]}
        ]
        lr_pairs = pd.DataFrame(
            [
                {
                    "pair_id": "lr1",
                    "ligand": "L",
                    "receptor": "R",
                    "programs": "exhaustion",
                    "pair_complete": False,
                }
            ]
        )
        result = find_program_overlap_combinations(sender_receiver_pairs, lr_pairs)
        assert result == []
