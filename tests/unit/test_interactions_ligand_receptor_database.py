"""Unit tests for xenium_tcr_ecology.interactions.ligand_receptor_database (`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.interactions.ligand_receptor_database import (
    CANDIDATE_LR_PAIRS,
    check_pair_panel_support,
)


class TestCheckPairPanelSupport:
    def test_both_present_gives_complete_pair(self):
        pair = {"pair_id": "X", "ligand": "GENE_A", "receptor": "GENE_B", "programs": []}
        result = check_pair_panel_support(pair, {"GENE_A", "GENE_B"})
        assert result["ligand_present"] is True
        assert result["receptor_present"] is True
        assert result["pair_complete"] is True

    def test_only_ligand_present_gives_incomplete_pair(self):
        pair = {"pair_id": "X", "ligand": "GENE_A", "receptor": "GENE_B", "programs": []}
        result = check_pair_panel_support(pair, {"GENE_A"})
        assert result["ligand_present"] is True
        assert result["receptor_present"] is False
        assert result["pair_complete"] is False

    def test_neither_present_gives_incomplete_pair(self):
        pair = {"pair_id": "X", "ligand": "GENE_A", "receptor": "GENE_B", "programs": []}
        result = check_pair_panel_support(pair, set())
        assert result["pair_complete"] is False


class TestCandidateLrPairs:
    def test_every_real_pair_has_a_ligand_and_receptor_gene(self):
        for pair in CANDIDATE_LR_PAIRS:
            assert pair["ligand"]
            assert pair["receptor"]
            assert pair["ligand"] != pair["receptor"]

    def test_pair_ids_are_unique(self):
        ids = [p["pair_id"] for p in CANDIDATE_LR_PAIRS]
        assert len(ids) == len(set(ids))
