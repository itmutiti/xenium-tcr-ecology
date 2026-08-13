"""Unit tests for xenium_tcr_ecology.preprocess.feature_classification (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.preprocess.feature_classification import classify_feature


class TestClassifyFeature:
    def test_classifies_standard_cdr3_probe(self):
        assert classify_feature("230322_CASSLEQGTQYF_TRB") == "tcr_cdr3_probe"

    def test_classifies_cdr3_probe_with_batch_letter_suffix(self):
        """Regression test: an earlier, narrower pattern (used in Phase
        4.08) required the date prefix to be exactly 6 digits, missing the
        real "231004B" batch (6 digits + a single letter) -- 17 real probes,
        confirmed present in P12 and P28's panels."""
        assert classify_feature("231004B_CASRDSPSTDTQYF_TRB") == "tcr_cdr3_probe"

    def test_classifies_tra_and_trb_suffixes(self):
        assert classify_feature("240501_CAVRDSNYQLIW_TRA") == "tcr_cdr3_probe"
        assert classify_feature("240501_CASSPGQGDTQYF_TRB") == "tcr_cdr3_probe"

    def test_classifies_hpv_probe(self):
        assert classify_feature("HPV16_E6") == "hpv_probe"
        assert classify_feature("HPV16_L1") == "hpv_probe"

    def test_classifies_standard_biological_gene(self):
        assert classify_feature("EPCAM") == "biological_gene"
        assert classify_feature("CD3D") == "biological_gene"

    def test_classifies_negative_control_probe(self):
        assert classify_feature("NegControlProbe_00022") == "negative_control_probe"

    def test_classifies_negative_control_codeword(self):
        assert classify_feature("NegControlCodeword_0508") == "negative_control_codeword"

    def test_classifies_unassigned_codeword(self):
        assert classify_feature("UnassignedCodeword_0069") == "unassigned_codeword"

    def test_does_not_misclassify_a_gene_that_merely_contains_tr_substring(self):
        # Sanity check against over-eager pattern matching: a real gene name
        # containing "TR" should not be swept into the CDR3 bucket.
        assert classify_feature("TRAC") == "biological_gene"
