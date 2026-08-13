"""Unit tests for xenium_tcr_ecology.tcr.cell_calls (`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.tcr.cell_calls import call_cell_detections


class TestCallCellDetections:
    def test_detects_own_patient_probe(self):
        counts = np.array([[1, 0]])  # one cell, two probes
        probe_names = ["probeA", "probeB"]
        probe_intended_patient = {"probeA": "P1", "probeB": "P2"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(counts, probe_names, probe_intended_patient, cell_patient_ids)
        assert result.loc[0, "n_probes_detected"] == 1
        assert result.loc[0, "detected_probes"] == "probeA"

    def test_ignores_detection_for_a_probe_not_intended_for_this_patient(self):
        # Cell is from P2 but has nonzero counts for probeA, which is
        # intended for P1 -- must not be called, even though the probe is
        # physically present (`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`, `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s batch-sharing finding).
        counts = np.array([[5, 0]])
        probe_names = ["probeA", "probeB"]
        probe_intended_patient = {"probeA": "P1", "probeB": "P2"}
        cell_patient_ids = pd.Series(["P2"])
        result = call_cell_detections(counts, probe_names, probe_intended_patient, cell_patient_ids)
        assert result.loc[0, "n_probes_detected"] == 0
        assert not result.loc[0, "any_detection"]

    def test_multi_probe_detection_is_flagged_ambiguous(self):
        counts = np.array([[1, 1]])
        probe_names = ["probeA", "probeC"]
        probe_intended_patient = {"probeA": "P1", "probeC": "P1"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(counts, probe_names, probe_intended_patient, cell_patient_ids)
        assert result.loc[0, "n_probes_detected"] == 2
        assert result.loc[0, "is_multi_probe_ambiguous"]

    def test_single_detection_is_not_ambiguous(self):
        counts = np.array([[1, 0]])
        probe_names = ["probeA", "probeC"]
        probe_intended_patient = {"probeA": "P1", "probeC": "P1"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(counts, probe_names, probe_intended_patient, cell_patient_ids)
        assert not result.loc[0, "is_multi_probe_ambiguous"]

    def test_n_probes_evaluated_counts_only_in_scope_probes(self):
        counts = np.array([[0, 0, 0]])
        probe_names = ["probeA", "probeB", "probeC"]
        probe_intended_patient = {"probeA": "P1", "probeB": "P2", "probeC": "P1"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(counts, probe_names, probe_intended_patient, cell_patient_ids)
        assert result.loc[0, "n_probes_evaluated"] == 2

    def test_zero_detection_is_a_real_negative_call(self):
        counts = np.array([[0, 0]])
        probe_names = ["probeA", "probeB"]
        probe_intended_patient = {"probeA": "P1", "probeB": "P1"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(counts, probe_names, probe_intended_patient, cell_patient_ids)
        assert not result.loc[0, "any_detection"]
        assert result.loc[0, "detected_probes"] == ""

    def test_tra_trb_pair_is_flagged_as_likely_single_clone(self):
        counts = np.array([[1, 1]])
        probe_names = ["probeA_TRA", "probeB_TRB"]
        probe_intended_patient = {"probeA_TRA": "P1", "probeB_TRB": "P1"}
        probe_chain = {"probeA_TRA": "TRA", "probeB_TRB": "TRB"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(
            counts, probe_names, probe_intended_patient, cell_patient_ids, probe_chain
        )
        assert result.loc[0, "likely_single_clone_tra_trb_pair"]
        assert not result.loc[0, "is_multi_probe_ambiguous_excluding_likely_pairs"]

    def test_two_tra_probes_is_not_a_likely_pair(self):
        counts = np.array([[1, 1]])
        probe_names = ["probeA_TRA", "probeB_TRA"]
        probe_intended_patient = {"probeA_TRA": "P1", "probeB_TRA": "P1"}
        probe_chain = {"probeA_TRA": "TRA", "probeB_TRA": "TRA"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(
            counts, probe_names, probe_intended_patient, cell_patient_ids, probe_chain
        )
        assert not result.loc[0, "likely_single_clone_tra_trb_pair"]
        assert result.loc[0, "is_multi_probe_ambiguous_excluding_likely_pairs"]

    def test_three_probes_is_not_a_likely_pair_even_with_tra_and_trb_present(self):
        counts = np.array([[1, 1, 1]])
        probe_names = ["probeA_TRA", "probeB_TRB", "probeC_TRA"]
        probe_intended_patient = {"probeA_TRA": "P1", "probeB_TRB": "P1", "probeC_TRA": "P1"}
        probe_chain = {"probeA_TRA": "TRA", "probeB_TRB": "TRB", "probeC_TRA": "TRA"}
        cell_patient_ids = pd.Series(["P1"])
        result = call_cell_detections(
            counts, probe_names, probe_intended_patient, cell_patient_ids, probe_chain
        )
        assert not result.loc[0, "likely_single_clone_tra_trb_pair"]
        assert result.loc[0, "is_multi_probe_ambiguous_excluding_likely_pairs"]
