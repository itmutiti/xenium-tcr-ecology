"""Unit tests for xenium_tcr_ecology.interactions.sender_receiver_pairs (`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.interactions.sender_receiver_pairs import (
    SENDER_RECEIVER_PAIRS,
    TGFB_GENE_SET,
    build_program_gene_sets_for_interactions,
    validate_sender_receiver_pairs,
)


class TestBuildProgramGeneSetsForInteractions:
    def test_includes_new_chemokine_and_tgf_beta_sets(self):
        result = build_program_gene_sets_for_interactions()
        assert "chemokine" in result
        assert "tgf_beta" in result
        assert len(result["chemokine"]) > 0

    def test_tgf_beta_set_is_really_empty_not_fabricated(self):
        assert TGFB_GENE_SET == []

    def test_includes_all_original_preprocessing_programs(self):
        result = build_program_gene_sets_for_interactions()
        for program in [
            "cytotoxicity",
            "exhaustion",
            "activation",
            "interferon",
            "proliferation",
            "antigen_presentation",
        ]:
            assert program in result


class TestValidateSenderReceiverPairs:
    def test_real_pairs_pass_validation_with_correct_lineages_and_programs(self):
        valid_lineages = {"Epithelial_Tumour", "Fibroblast", "Myeloid", "Dendritic_cell", "T_cell"}
        valid_programs = set(build_program_gene_sets_for_interactions().keys())
        errors = validate_sender_receiver_pairs(
            SENDER_RECEIVER_PAIRS, valid_lineages, valid_programs
        )
        assert errors == []

    def test_unknown_sender_lineage_is_caught(self):
        pairs = [
            {
                "pair_id": "bad",
                "sender": "NotALineage",
                "receiver": "T_cell",
                "relevant_programs": [],
            }
        ]
        errors = validate_sender_receiver_pairs(pairs, {"T_cell"}, set())
        assert len(errors) == 1
        assert "NotALineage" in errors[0]

    def test_unknown_program_is_caught(self):
        pairs = [
            {
                "pair_id": "bad",
                "sender": "T_cell",
                "receiver": "T_cell",
                "relevant_programs": ["not_a_real_program"],
            }
        ]
        errors = validate_sender_receiver_pairs(pairs, {"T_cell"}, {"exhaustion"})
        assert len(errors) == 1
        assert "not_a_real_program" in errors[0]
