"""Unit tests for xenium_tcr_ecology.tcr.vdj_ground_truth_validation (`08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_vdj_ground_truth.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.tcr.vdj_ground_truth_validation import (
    assign_cell_hashes,
    load_aggregation_order,
    load_pool_to_patients,
    pool_key_from_sample_name,
    render_ground_truth_validation_report,
)


class TestLoadAggregationOrder:
    def test_real_suffix_to_sample_mapping_is_one_indexed_by_row_order(self, tmp_path):
        path = tmp_path / "aggregation.csv"
        path.write_text("sample_id,molecule_h5\nchip1pool3,x\nchip2pool2,y\nchip1pool1,z\n")
        result = load_aggregation_order(path)
        assert result == {1: "chip1pool3", 2: "chip2pool2", 3: "chip1pool1"}


class TestLoadPoolToPatients:
    def test_real_pool_to_hash_to_patient_mapping(self, tmp_path):
        path = tmp_path / "patient_matrix.txt"
        path.write_text(
            '""\t"pool1"\t"pool2"\n"hash1"\t1\t6\n"hash2"\t8\t9\n"hash3"\t15\t16\n"hash4"\t28\t24\n'
        )
        result = load_pool_to_patients(path)
        assert result["pool1"] == {"hash1": 1, "hash2": 8, "hash3": 15, "hash4": 28}
        assert result["pool2"] == {"hash1": 6, "hash2": 9, "hash3": 16, "hash4": 24}


class TestPoolKeyFromSampleName:
    def test_real_chip1_prefix_stripped(self):
        assert pool_key_from_sample_name("chip1pool3") == "pool3"

    def test_real_chip2_prefix_stripped(self):
        assert pool_key_from_sample_name("chip2pool16") == "pool16"

    def test_real_name_without_chip_prefix_unchanged(self):
        assert pool_key_from_sample_name("pool5") == "pool5"


class TestAssignCellHashes:
    def test_real_confident_top_hash_assigned(self):
        counts = np.array([[100, 5, 2, 1]])
        result = assign_cell_hashes(counts)
        assert result[0] == 0

    def test_real_ambiguous_cell_gets_unassigned(self):
        # top and second-highest too close -- below MIN_HASH_CONFIDENCE_RATIO
        counts = np.array([[50, 40, 2, 1]])
        result = assign_cell_hashes(counts)
        assert result[0] == -1

    def test_real_low_absolute_count_gets_unassigned_even_if_dominant(self):
        # top count below MIN_HASH_UMI_COUNT despite a clean ratio
        counts = np.array([[5, 0, 0, 0]])
        result = assign_cell_hashes(counts)
        assert result[0] == -1

    def test_real_mixed_batch_assigns_each_cell_independently(self):
        counts = np.array([[100, 5, 2, 1], [1, 200, 3, 2], [50, 40, 2, 1]])
        result = assign_cell_hashes(counts)
        assert list(result) == [0, 1, -1]


class TestRenderGroundTruthValidationReport:
    def test_real_report_written(self, tmp_path):
        comparison = pd.DataFrame(
            {
                "probe_name": ["p1", "p2", "p3", "p4"],
                "intended_patient": ["P01", "P01", "P09", "P09"],
                "tcr_chain": ["TRA", "TRB", "TRA", "TRB"],
                "cdr3_amino_acid_sequence": ["CAA", "CAB", "CAC", "CAD"],
                "xenium_detection_rate": [0.01, 0.02, 0.005, 0.03],
                "found_in_real_vdj_ground_truth": [True, True, True, False],
                "vdj_ground_truth_n_cells": [100, 50, 20, 0],
                "vdj_ground_truth_rank_within_patient_chain": [1.0, 2.0, 3.0, float("nan")],
            }
        )
        output_path = tmp_path / "reports" / "tcr" / "probe_vdj_ground_truth_validation.pdf"
        render_ground_truth_validation_report(comparison, output_path)
        assert output_path.is_file()
        assert output_path.stat().st_size > 0
