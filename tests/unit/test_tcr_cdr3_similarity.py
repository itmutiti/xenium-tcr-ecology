"""Unit tests for xenium_tcr_ecology.tcr.cdr3_similarity (`08_tcr_clonal_analysis/05_screen_cdr3_cross_patient_similarity.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.tcr.cdr3_similarity import (
    levenshtein_distance,
    screen_cdr3_pairwise_similarity,
)


class TestLevenshteinDistance:
    def test_identical_strings_have_distance_zero(self):
        assert levenshtein_distance("CASSLEG", "CASSLEG") == 0

    def test_single_substitution_has_distance_one(self):
        assert levenshtein_distance("CASSLEG", "CASSLEA") == 1

    def test_single_insertion_has_distance_one(self):
        assert levenshtein_distance("CASSLEG", "CASSLEGA") == 1

    def test_single_deletion_has_distance_one(self):
        assert levenshtein_distance("CASSLEG", "CASSLE") == 1

    def test_empty_string_distance_is_length_of_other(self):
        assert levenshtein_distance("", "ABCDE") == 5
        assert levenshtein_distance("ABCDE", "") == 5

    def test_completely_different_strings(self):
        assert levenshtein_distance("AAAA", "GGGG") == 4


class TestScreenCdr3PairwiseSimilarity:
    def _registry(self):
        return pd.DataFrame(
            {
                "probe_name": ["p1", "p2", "p3", "p4"],
                "cdr3_amino_acid_sequence": [
                    "CASSLEGATDTQYF",
                    "CASSLEGATDTQYA",
                    "CAAQNSGYSTLTF",
                    "CASSLEGATDT",
                ],
                "tcr_chain": ["TRB", "TRB", "TRA", "TRB"],
                "patients_with_probe": ["P1", "P2", "P1;P2", "P3"],
            }
        )

    def test_flags_near_identical_same_chain_pair(self):
        result = screen_cdr3_pairwise_similarity(self._registry(), max_distance=2)
        pair = result[(result["probe_a"] == "p1") & (result["probe_b"] == "p2")]
        assert len(pair) == 1
        assert pair.iloc[0]["edit_distance"] == 1

    def test_does_not_compare_different_chains(self):
        result = screen_cdr3_pairwise_similarity(self._registry(), max_distance=10)
        cross_chain = result[((result["probe_a"] == "p3") | (result["probe_b"] == "p3"))]
        assert len(cross_chain) == 0

    def test_cross_patient_flag_is_correct(self):
        result = screen_cdr3_pairwise_similarity(self._registry(), max_distance=2)
        pair = result[(result["probe_a"] == "p1") & (result["probe_b"] == "p2")]
        assert pair.iloc[0]["is_cross_patient"]  # P1 vs P2, disjoint

    def test_shared_patient_is_not_cross_patient(self):
        registry = pd.DataFrame(
            {
                "probe_name": ["p1", "p2"],
                "cdr3_amino_acid_sequence": ["CASSLEGATDTQYF", "CASSLEGATDTQYA"],
                "tcr_chain": ["TRB", "TRB"],
                "patients_with_probe": ["P1;P2", "P2;P3"],
            }
        )
        result = screen_cdr3_pairwise_similarity(registry, max_distance=2)
        assert not result.iloc[0]["is_cross_patient"]  # both include P2

    def test_dissimilar_pairs_are_not_flagged(self):
        result = screen_cdr3_pairwise_similarity(self._registry(), max_distance=2)
        pair = result[(result["probe_a"] == "p1") & (result["probe_b"] == "p4")]
        assert len(pair) == 0  # distance > 2
