"""Unit tests for xenium_tcr_ecology.interactions.barrier_pathways (`14_spatial_interactions_and_barriers/04_analyse_barrier_pathways.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.interactions.barrier_pathways import (
    classify_clone_engagement,
    compare_interface_programs,
    identify_barrier_interface_cells,
)


class TestClassifyCloneEngagement:
    def test_real_median_split_below_is_excluded(self):
        engagement = pd.Series([0.1, 0.2, 0.8, 0.9])
        result = classify_clone_engagement(engagement)
        assert list(result) == ["excluded", "excluded", "engaged", "engaged"]

    def test_real_value_exactly_at_median_is_engaged_not_excluded(self):
        # Median split convention: values >= median are "engaged" --
        # checked explicitly so a boundary value is never silently
        # misclassified.
        engagement = pd.Series([0.1, 0.5, 0.5, 0.9])
        result = classify_clone_engagement(engagement)
        assert list(result) == ["excluded", "engaged", "engaged", "engaged"]


class TestIdentifyBarrierInterfaceCells:
    def test_real_interface_cells_restricted_to_target_barrier_groups(self):
        # Path: clone cell (idx 3) -> fibroblast (idx 2) -> vascular (idx 1) -> tumour (idx 0).
        predecessors = np.array([-9999, 0, 1, 2])
        barrier_group = np.array(["tumour", "vascular", "fibroblast", "other"], dtype=object)
        result = identify_barrier_interface_cells(
            predecessors,
            tumour_idx_set={0},
            clone_cell_indices=[3],
            barrier_group=barrier_group,
            target_groups=("fibroblast", "suppressive_myeloid"),
        )
        assert result == {2}

    def test_real_no_interface_cells_when_path_is_empty(self):
        predecessors = np.array([-9999, 0])
        barrier_group = np.array(["tumour", "other"], dtype=object)
        result = identify_barrier_interface_cells(
            predecessors,
            tumour_idx_set={0},
            clone_cell_indices=[1],
            barrier_group=barrier_group,
            target_groups=("fibroblast", "suppressive_myeloid"),
        )
        assert result == set()

    def test_real_multiple_clone_cells_union_their_interface_cells(self):
        # Clone cell A (idx 2) -> fibroblast (idx 1) -> tumour (idx 0).
        # Clone cell B (idx 4) -> suppressive_myeloid (idx 3) -> tumour (idx 0).
        predecessors = np.array([-9999, 0, 1, 0, 3])
        barrier_group = np.array(
            ["tumour", "fibroblast", "other", "suppressive_myeloid", "other"], dtype=object
        )
        result = identify_barrier_interface_cells(
            predecessors,
            tumour_idx_set={0},
            clone_cell_indices=[2, 4],
            barrier_group=barrier_group,
            target_groups=("fibroblast", "suppressive_myeloid"),
        )
        assert result == {1, 3}


class TestCompareInterfacePrograms:
    def test_real_elevated_excluded_group_score_is_detected(self):
        rng = np.random.default_rng(0)
        n = 60
        interface_table = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(n)],
                "clone_class": ["excluded"] * (n // 2) + ["engaged"] * (n // 2),
            }
        )
        program_scores = pd.DataFrame(
            {
                "checkpoint_score": np.concatenate(
                    [rng.normal(2.0, 0.5, n // 2), rng.normal(0.0, 0.5, n // 2)]
                ),
            },
            index=interface_table["cell_id"],
        )
        result = compare_interface_programs(interface_table, program_scores, ["checkpoint_score"])
        row = result.iloc[0]
        assert row["mean_excluded"] > row["mean_engaged"]
        assert row["pvalue"] < 0.05

    def test_real_no_signal_gives_non_significant_pvalue(self):
        rng = np.random.default_rng(1)
        n = 60
        interface_table = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(n)],
                "clone_class": ["excluded"] * (n // 2) + ["engaged"] * (n // 2),
            }
        )
        program_scores = pd.DataFrame(
            {"checkpoint_score": rng.normal(0.0, 1.0, n)},
            index=interface_table["cell_id"],
        )
        result = compare_interface_programs(interface_table, program_scores, ["checkpoint_score"])
        assert result.iloc[0]["pvalue"] > 0.05

    def test_bh_adjustment_is_never_smaller_than_raw_pvalue(self):
        rng = np.random.default_rng(2)
        n = 40
        interface_table = pd.DataFrame(
            {
                "cell_id": [f"c{i}" for i in range(n)],
                "clone_class": ["excluded"] * (n // 2) + ["engaged"] * (n // 2),
            }
        )
        program_scores = pd.DataFrame(
            {
                "a_score": rng.normal(0, 1, n),
                "b_score": rng.normal(0, 1, n),
            },
            index=interface_table["cell_id"],
        )
        result = compare_interface_programs(interface_table, program_scores, ["a_score", "b_score"])
        assert (result["pvalue_bh"] >= result["pvalue"] - 1e-12).all()
