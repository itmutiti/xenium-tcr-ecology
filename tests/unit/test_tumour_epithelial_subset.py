"""Unit tests for xenium_tcr_ecology.tumour.epithelial_subset (`07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.epithelial_subset import (
    EPITHELIAL_LINEAGE_LABEL,
    compute_cell_weighted_patient_dominance,
    compute_joint_cluster_patient_dominance,
    subset_epithelial_cells,
)


def _make_adata(cell_ids, patient_ids):
    n = len(cell_ids)
    return ad.AnnData(
        X=np.zeros((n, 3), dtype=np.float32),
        obs=pd.DataFrame({"patient_id": patient_ids}, index=cell_ids),
    )


class TestSubsetEpithelialCells:
    def test_subsets_only_epithelial_lineage(self):
        cell_ids = ["c1", "c2", "c3"]
        adata = _make_adata(cell_ids, ["P1", "P1", "P2"])
        final_annotations = pd.DataFrame(
            {
                "final_lineage": [EPITHELIAL_LINEAGE_LABEL, "T_cell", EPITHELIAL_LINEAGE_LABEL],
                "confidence": [0.9, 0.5, 0.3],
                "is_ambiguous": [False, False, True],
            },
            index=cell_ids,
        )
        sub = subset_epithelial_cells(adata, final_annotations)
        assert set(sub.obs_names) == {"c1", "c3"}
        assert sub.obs.loc["c3", "is_ambiguous"]

    def test_raises_when_no_epithelial_cells_in_annotations(self):
        adata = _make_adata(["c1"], ["P1"])
        final_annotations = pd.DataFrame(
            {"final_lineage": ["T_cell"], "confidence": [0.5], "is_ambiguous": [False]},
            index=["c1"],
        )
        with pytest.raises(PipelineError, match=EPITHELIAL_LINEAGE_LABEL):
            subset_epithelial_cells(adata, final_annotations)

    def test_raises_when_no_overlap_with_adata(self):
        adata = _make_adata(["c1"], ["P1"])
        final_annotations = pd.DataFrame(
            {
                "final_lineage": [EPITHELIAL_LINEAGE_LABEL],
                "confidence": [0.5],
                "is_ambiguous": [False],
            },
            index=["other_cell"],
        )
        with pytest.raises(PipelineError, match="overlap"):
            subset_epithelial_cells(adata, final_annotations)


class TestComputeJointClusterPatientDominance:
    def test_fully_patient_dominated_clusters_score_one(self):
        clusters = pd.Series(["0", "0", "1", "1"])
        patients = pd.Series(["P1", "P1", "P2", "P2"])
        assert compute_joint_cluster_patient_dominance(clusters, patients) == pytest.approx(1.0)

    def test_evenly_mixed_cluster_scores_low(self):
        clusters = pd.Series(["0", "0", "0", "0"])
        patients = pd.Series(["P1", "P2", "P3", "P4"])
        assert compute_joint_cluster_patient_dominance(clusters, patients) == pytest.approx(0.25)

    def test_partial_dominance(self):
        # Cluster '0': 3xP1, 1xP2 -> dominance 0.75; cluster '1': all P3 -> 1.0.
        clusters = pd.Series(["0", "0", "0", "0", "1", "1"])
        patients = pd.Series(["P1", "P1", "P1", "P2", "P3", "P3"])
        assert compute_joint_cluster_patient_dominance(clusters, patients) == pytest.approx(
            (0.75 + 1.0) / 2
        )


class TestComputeCellWeightedPatientDominance:
    def test_equals_unweighted_when_clusters_are_equal_size(self):
        clusters = pd.Series(["0", "0", "1", "1"])
        patients = pd.Series(["P1", "P1", "P2", "P2"])
        unweighted = compute_joint_cluster_patient_dominance(clusters, patients)
        weighted = compute_cell_weighted_patient_dominance(clusters, patients)
        assert weighted == pytest.approx(unweighted)

    def test_large_dominant_cluster_pulls_weighted_average_up(self):
        # A large, fully patient-exclusive cluster ('0', 100 cells) plus a
        # small, evenly-mixed cluster ('1', 4 cells): the cluster-unweighted
        # mean treats both clusters equally ((1.0 + 0.25) / 2 = 0.625), but
        # most CELLS are in the highly patient-dominant cluster, so the
        # cell-weighted mean should be much closer to 1.0.
        clusters = pd.Series(["0"] * 100 + ["1"] * 4)
        patients = pd.Series(["P1"] * 100 + ["P2", "P3", "P4", "P5"])
        unweighted = compute_joint_cluster_patient_dominance(clusters, patients)
        weighted = compute_cell_weighted_patient_dominance(clusters, patients)
        assert unweighted == pytest.approx(0.625)
        assert weighted > unweighted
        assert weighted == pytest.approx((100 * 1.0 + 4 * 0.25) / 104)
