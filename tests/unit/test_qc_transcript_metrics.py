"""Unit tests for xenium_tcr_ecology.qc.transcript_metrics (`04_quality_control/01_compute_transcript_level_qc_metrics.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.qc.transcript_metrics import _classify_feature_name, _is_assigned_to_cell


class TestIsAssignedToCell:
    def test_real_cell_id_is_assigned(self):
        assert _is_assigned_to_cell("aaadggoi-1") is True

    def test_literal_unassigned_sentinel_is_not_assigned(self):
        """Regression test: real Xenium data uses the literal string
        "UNASSIGNED" (not an empty string) to mark a transcript with no
        assigned cell -- confirmed against real data (27.4% of one section's
        transcripts). An earlier version of this module checked
        `cell_id != ""` only, which silently treated every "UNASSIGNED"
        transcript as assigned, producing n_unassigned == 0 for all 18
        sections in the resulting report."""
        assert _is_assigned_to_cell("UNASSIGNED") is False

    def test_empty_string_is_not_assigned(self):
        assert _is_assigned_to_cell("") is False


class TestClassifyFeatureName:
    def test_classifies_negative_control_probe(self):
        assert _classify_feature_name("NegControlProbe_00022") == "negative_control_probe"

    def test_classifies_negative_control_codeword(self):
        assert _classify_feature_name("NegControlCodeword_0508") == "negative_control_codeword"

    def test_classifies_unassigned_codeword(self):
        assert _classify_feature_name("UnassignedCodeword_0069") == "unassigned_codeword"

    def test_classifies_gene_and_custom_probes_as_gene_expression(self):
        assert _classify_feature_name("EPCAM") == "gene_expression"
        assert _classify_feature_name("HPV16_E6") == "gene_expression"
        assert _classify_feature_name("230322_CASSLEQGTQYF_TRB") == "gene_expression"
