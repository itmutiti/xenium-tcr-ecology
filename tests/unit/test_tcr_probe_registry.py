"""Unit tests for xenium_tcr_ecology.tcr.probe_registry (`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.tcr.probe_registry import parse_cdr3_probe_name


class TestParseCdr3ProbeName:
    def test_parses_probe_without_batch_letter(self):
        result = parse_cdr3_probe_name("230322_CAAQNSGYSTLTF_TRA")
        assert result == {
            "probe_name": "230322_CAAQNSGYSTLTF_TRA",
            "date_batch_prefix": "230322",
            "cdr3_amino_acid_sequence": "CAAQNSGYSTLTF",
            "tcr_chain": "TRA",
        }

    def test_parses_probe_with_batch_letter(self):
        result = parse_cdr3_probe_name("231004B_CASRDSPSTDTQYF_TRB")
        assert result["date_batch_prefix"] == "231004B"
        assert result["tcr_chain"] == "TRB"
        assert result["cdr3_amino_acid_sequence"] == "CASRDSPSTDTQYF"

    def test_non_matching_name_returns_none(self):
        assert parse_cdr3_probe_name("EPCAM") is None
        assert parse_cdr3_probe_name("HPV16_E6") is None

    def test_malformed_chain_suffix_returns_none(self):
        assert parse_cdr3_probe_name("230322_CAAQNSGYSTLTF_TRG") is None
