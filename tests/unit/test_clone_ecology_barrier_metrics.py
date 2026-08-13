"""Unit tests for xenium_tcr_ecology.clone_ecology.barrier_metrics (`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.clone_ecology.barrier_metrics import (
    assign_barrier_group,
    compute_barrier_fraction,
    compute_empirical_pvalue_greater,
    trace_intermediate_path,
)


class TestAssignBarrierGroup:
    def test_fibroblast_lineage_maps_to_fibroblast_group(self):
        lineage = np.array(["Fibroblast", "T_cell"])
        substate = np.array(["Activated_CAF", None])
        result = assign_barrier_group(lineage, substate)
        assert list(result) == ["fibroblast", "other"]

    def test_vascular_lineages_map_to_vascular_group(self):
        lineage = np.array(["Endothelial", "Perivascular_SmoothMuscle", "T_cell"])
        substate = np.array([None, None, None])
        result = assign_barrier_group(lineage, substate)
        assert list(result) == ["vascular", "vascular", "other"]

    def test_macrophage_substate_maps_to_suppressive_myeloid(self):
        lineage = np.array(["Myeloid", "Myeloid"])
        substate = np.array(["Macrophage", "Monocyte"])
        result = assign_barrier_group(lineage, substate)
        assert list(result) == ["suppressive_myeloid", "other"]


class TestTraceIntermediatePath:
    def test_direct_neighbour_has_no_intermediate_cells(self):
        # predecessors[3] = 2 (a tumour cell) -- cell 3 is directly adjacent to tumour.
        predecessors = np.array([-9999, -9999, -9999, 2])
        result = trace_intermediate_path(predecessors, target_idx=3, source_idx_set={2})
        assert result == []

    def test_two_hop_path_has_one_intermediate_cell(self):
        # 3 -> 1 (intermediate) -> 0 (tumour source)
        predecessors = np.array([-9999, 0, -9999, 1])
        result = trace_intermediate_path(predecessors, target_idx=3, source_idx_set={0})
        assert result == [1]

    def test_unreachable_target_returns_empty(self):
        predecessors = np.array([-9999, -9999])
        result = trace_intermediate_path(predecessors, target_idx=1, source_idx_set={0})
        assert result == []


class TestComputeBarrierFraction:
    def test_fraction_of_matching_group(self):
        groups = np.array(["fibroblast", "fibroblast", "vascular", "other"])
        assert compute_barrier_fraction(groups, "fibroblast") == 0.5

    def test_empty_path_gives_nan(self):
        assert np.isnan(compute_barrier_fraction(np.array([]), "fibroblast"))


class TestComputeEmpiricalPvalueGreater:
    def test_observed_far_above_null_gives_small_pvalue(self):
        null_values = np.zeros(199)
        result = compute_empirical_pvalue_greater(observed=1.0, null_values=null_values)
        assert result == 1.0 / 200.0

    def test_observed_far_below_null_gives_pvalue_near_one(self):
        null_values = np.ones(199)
        result = compute_empirical_pvalue_greater(observed=0.0, null_values=null_values)
        assert result == 1.0
