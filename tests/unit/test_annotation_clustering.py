"""Unit tests for xenium_tcr_ecology.annotation.clustering (`06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.clustering import (
    compute_joint_clustering,
    compute_within_patient_clustering,
    select_clustering_gene_mask,
)


class TestSelectClusteringGeneMask:
    def test_selects_only_biological_genes(self):
        var_names = pd.Index(["GENE1", "GENE2", "230322_CASSLEQGTQYF_TRB"])
        feature_annotation = pd.DataFrame(
            {
                "feature_name": var_names,
                "feature_class": ["biological_gene", "biological_gene", "tcr_cdr3_probe"],
            }
        )
        mask = select_clustering_gene_mask(var_names, feature_annotation)
        assert mask.tolist() == [True, True, False]

    def test_raises_if_no_biological_genes(self):
        var_names = pd.Index(["A"])
        feature_annotation = pd.DataFrame(
            {"feature_name": ["A"], "feature_class": ["tcr_cdr3_probe"]}
        )
        with pytest.raises(PipelineError, match="No biological_gene"):
            select_clustering_gene_mask(var_names, feature_annotation)


def _make_adata(n_cells_per_group=80, n_genes=20, rng_seed=0):
    """Two well-separated synthetic expression 'blobs' plus a patient
    label, so clustering has an obvious, checkable answer."""
    rng = np.random.default_rng(rng_seed)
    n_patients = 3

    group_a = rng.poisson(2, size=(n_cells_per_group, n_genes)).astype(np.float32)
    group_b = rng.poisson(2, size=(n_cells_per_group, n_genes)).astype(np.float32)
    group_b[:, :5] += 15  # first 5 genes strongly elevated -> a distinct cluster
    X = np.vstack([group_a, group_b])

    adata = ad.AnnData(X=X)
    adata.obs_names = [f"cell{i}" for i in range(X.shape[0])]
    adata.var_names = [f"GENE{i}" for i in range(n_genes)]
    adata.obs["patient_id"] = [f"P{(i % n_patients) + 1}" for i in range(X.shape[0])]
    adata.layers["lognorm"] = np.log1p(X)
    return adata


class TestComputeJointClustering:
    def test_recovers_the_two_synthetic_groups(self):
        adata = _make_adata()
        gene_mask = np.ones(adata.n_vars, dtype=bool)
        result = compute_joint_clustering(
            adata, layer="lognorm", gene_mask=gene_mask, resolutions=[0.5], n_top_genes=20
        )

        col = "joint_leiden_res0.5"
        assert col in result.columns
        # The two synthetic groups (first half vs second half of cells)
        # should be predominantly assigned to different clusters.
        first_half_mode = result[col].iloc[:80].mode()[0]
        second_half_mode = result[col].iloc[80:].mode()[0]
        assert first_half_mode != second_half_mode

    def test_multiple_resolutions_produce_a_column_each(self):
        adata = _make_adata()
        gene_mask = np.ones(adata.n_vars, dtype=bool)
        result = compute_joint_clustering(
            adata, layer="lognorm", gene_mask=gene_mask, resolutions=[0.3, 1.0], n_top_genes=20
        )
        assert "joint_leiden_res0.3" in result.columns
        assert "joint_leiden_res1.0" in result.columns
        assert len(result) == adata.n_obs


class TestComputeWithinPatientClustering:
    def test_returns_a_label_per_cell_for_patients_above_min_size(self):
        adata = _make_adata()
        gene_mask = np.ones(adata.n_vars, dtype=bool)
        result = compute_within_patient_clustering(
            adata, layer="lognorm", gene_mask=gene_mask, min_cells=10, n_top_genes=20
        )
        assert result.notna().all()
        # Labels are prefixed with patient_id to keep them globally unique.
        assert all(str(v).startswith(pid) for v, pid in zip(result, adata.obs["patient_id"]))

    def test_skips_patients_below_min_cells(self):
        adata = _make_adata()
        gene_mask = np.ones(adata.n_vars, dtype=bool)
        result = compute_within_patient_clustering(
            adata, layer="lognorm", gene_mask=gene_mask, min_cells=10_000, n_top_genes=20
        )
        assert result.isna().all()
