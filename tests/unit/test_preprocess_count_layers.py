"""Unit tests for xenium_tcr_ecology.preprocess.count_layers (`05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.count_layers import (
    build_analysis_count_layers,
    compute_pearson_residuals,
    compute_size_factors,
    normalize_log1p,
)


def _make_feature_annotation(var_names, bio_gene_names):
    return pd.DataFrame(
        {
            "feature_name": var_names,
            "feature_class": [
                "biological_gene" if g in bio_gene_names else "tcr_cdr3_probe" for g in var_names
            ],
            "in_analysis_matrix": True,
        }
    )


class TestComputeSizeFactors:
    def test_size_factor_uses_only_exposure_genes(self):
        # 2 cells, 3 genes: gene 2 (index 2) is excluded from exposure and
        # has wildly different counts -- must not affect the size factor.
        X = np.array([[10.0, 10.0, 1000.0], [10.0, 10.0, 0.0]])
        exposure_mask = np.array([True, True, False])
        size_factors, target_sum = compute_size_factors(X, exposure_mask)
        # Both cells have identical exposure-gene totals (20), so their size
        # factors must be identical despite the huge difference in gene 2.
        assert size_factors[0] == size_factors[1]

    def test_target_sum_is_median_of_exposure_totals(self):
        X = np.array([[10.0], [20.0], [30.0]])
        exposure_mask = np.array([True])
        _, target_sum = compute_size_factors(X, exposure_mask)
        assert target_sum == 20.0

    def test_raises_if_all_exposure_counts_zero(self):
        X = np.zeros((3, 2))
        exposure_mask = np.array([True, False])
        with pytest.raises(PipelineError, match="zero counts"):
            compute_size_factors(X, exposure_mask)


class TestNormalizeLog1p:
    def test_dense_and_sparse_agree(self):
        X_dense = np.array([[4.0, 0.0], [0.0, 8.0]])
        X_sparse = sparse.csr_matrix(X_dense)
        size_factors = np.array([2.0, 4.0])

        dense_result = normalize_log1p(X_dense, size_factors)
        sparse_result = normalize_log1p(X_sparse, size_factors)

        np.testing.assert_allclose(dense_result, np.asarray(sparse_result.todense()), rtol=1e-6)

    def test_known_value(self):
        X = np.array([[10.0]])
        size_factors = np.array([2.0])
        result = normalize_log1p(X, size_factors)
        assert result[0, 0] == pytest.approx(np.log1p(5.0))


class TestComputePearsonResiduals:
    def test_matches_scanpy_when_exposure_is_every_gene(self):
        """Cross-validation: with the exposure mask set to all genes, this
        module's generalised formula must exactly reproduce scanpy's own
        sc.experimental.pp.normalize_pearson_residuals -- the strongest
        available correctness check on the manual reimplementation."""
        rng = np.random.default_rng(0)
        X = rng.poisson(3, size=(30, 8)).astype(np.float32)

        adata = ad.AnnData(X=X.copy())
        sc.experimental.pp.normalize_pearson_residuals(
            adata, theta=100, clip=None, check_values=False
        )
        expected = adata.X

        exposure_mask = np.ones(8, dtype=bool)
        result = compute_pearson_residuals(X, exposure_mask, theta=100, clip=None)

        np.testing.assert_allclose(result, expected, rtol=1e-5, atol=1e-5)

    def test_excluding_a_dominant_gene_changes_other_genes_residuals(self):
        """A cell with an extreme count in an excluded gene must not have
        its other genes' residuals distorted by that gene's exposure
        contribution -- the whole point of restricting the exposure basis."""
        X = np.array(
            [
                [10.0, 10.0, 10.0],
                [10.0, 10.0, 10.0],
                [10.0, 10.0, 5000.0],  # cell 2 has an extreme value in gene 2
            ]
        )
        exposure_mask = np.array([True, True, False])
        residuals_restricted = compute_pearson_residuals(X, exposure_mask, theta=100, clip=np.inf)

        full_mask = np.array([True, True, True])
        residuals_full = compute_pearson_residuals(X, full_mask, theta=100, clip=np.inf)

        # Cell 2's residual for gene 0 should be near-zero (unrestricted
        # exposure, matching its peers) when gene 2 is excluded from
        # exposure, but would be pulled negative if gene 2's huge count
        # were allowed to inflate cell 2's apparent size.
        assert abs(residuals_restricted[2, 0]) < abs(residuals_full[2, 0])


class TestBuildAnalysisCountLayers:
    def test_creates_all_four_layers(self):
        X = sparse.csr_matrix(np.array([[5.0, 0.0, 1.0], [0.0, 3.0, 0.0]], dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = ["GENE1", "GENE2", "230322_CASSLEQGTQYF_TRB"]
        feature_annotation = _make_feature_annotation(adata.var_names, {"GENE1", "GENE2"})

        result = build_analysis_count_layers(adata, feature_annotation)

        assert set(result.layers.keys()) == {"counts", "lognorm", "pearson_residuals", "detected"}
        np.testing.assert_array_equal(
            np.asarray(result.layers["counts"].todense()), np.asarray(X.todense())
        )
        np.testing.assert_array_equal(
            np.asarray(result.layers["detected"].todense()),
            np.array([[1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
        )
        assert result.var["is_exposure_gene"].tolist() == [True, True, False]

    def test_raises_if_no_biological_gene_features(self):
        X = sparse.csr_matrix(np.ones((2, 2), dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = ["A", "B"]
        feature_annotation = pd.DataFrame(
            {
                "feature_name": ["A", "B"],
                "feature_class": ["tcr_cdr3_probe", "tcr_cdr3_probe"],
                "in_analysis_matrix": True,
            }
        )
        with pytest.raises(PipelineError, match="biological_gene"):
            build_analysis_count_layers(adata, feature_annotation)

    def test_raises_on_missing_feature_annotation_column(self):
        X = sparse.csr_matrix(np.ones((2, 2), dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = ["A", "B"]
        with pytest.raises(PipelineError, match="missing required column"):
            build_analysis_count_layers(adata, pd.DataFrame({"feature_name": ["A", "B"]}))
