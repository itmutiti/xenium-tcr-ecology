"""Unit tests for xenium_tcr_ecology.tcr.resolve_calls (`08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.tcr.resolve_calls import classify_calls


def _row(any_detection, ambiguous, detected_probes):
    return {
        "any_detection": any_detection,
        "is_multi_probe_ambiguous_excluding_likely_pairs": ambiguous,
        "detected_probes": detected_probes,
    }


class TestClassifyCalls:
    def test_no_detection_is_unassigned(self):
        calls = pd.DataFrame([_row(False, False, "")])
        result = classify_calls(calls, probe_fpr={})
        assert result.iloc[0] == "unassigned"

    def test_genuinely_ambiguous_is_probable_multiplet(self):
        calls = pd.DataFrame([_row(True, True, "probeA;probeB")])
        result = classify_calls(calls, probe_fpr={"probeA": 0.1, "probeB": 0.1})
        assert result.iloc[0] == "probable_multiplet"

    def test_clean_singlet_with_reliable_probe(self):
        calls = pd.DataFrame([_row(True, False, "probeA")])
        result = classify_calls(calls, probe_fpr={"probeA": 0.1})
        assert result.iloc[0] == "singlet"

    def test_singlet_with_unreliable_probe_is_low_confidence(self):
        calls = pd.DataFrame([_row(True, False, "probeA")])
        result = classify_calls(calls, probe_fpr={"probeA": 0.9})
        assert result.iloc[0] == "low_confidence"

    def test_fpr_threshold_is_a_strict_greater_than(self):
        calls = pd.DataFrame([_row(True, False, "probeA")])
        result = classify_calls(calls, probe_fpr={"probeA": 0.5}, fpr_threshold=0.5)
        assert result.iloc[0] == "singlet"

    def test_missing_fpr_defaults_to_reliable(self):
        calls = pd.DataFrame([_row(True, False, "probeA")])
        result = classify_calls(calls, probe_fpr={})
        assert result.iloc[0] == "singlet"

    def test_multiplet_precedence_over_low_confidence(self):
        # Ambiguous AND involves an unreliable probe --
        # multiplet status takes precedence in the classification order.
        calls = pd.DataFrame([_row(True, True, "probeA;probeB")])
        result = classify_calls(calls, probe_fpr={"probeA": 0.9, "probeB": 0.9})
        assert result.iloc[0] == "probable_multiplet"

    def test_likely_pair_with_one_unreliable_probe_is_low_confidence(self):
        # Not flagged multiplet (excluded as a likely TRA/TRB pair
        # upstream, so is_multi_probe_ambiguous_excluding_likely_pairs is
        # False), but one of the two probes is unreliable.
        calls = pd.DataFrame([_row(True, False, "probeA_TRA;probeB_TRB")])
        result = classify_calls(calls, probe_fpr={"probeA_TRA": 0.1, "probeB_TRB": 0.9})
        assert result.iloc[0] == "low_confidence"
