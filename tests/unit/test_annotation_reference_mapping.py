"""Unit tests for xenium_tcr_ecology.annotation.reference_mapping (`06_cell_type_annotation/03_map_external_scrna_reference.py`)."""

from __future__ import annotations

import gzip

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.reference_mapping import (
    classify_by_nearest_centroid,
    compute_state_centroids,
    load_reference_matrix,
    validate_t_cell_state_markers,
    benchmark_transfer_confidence_degradation,
)


class TestValidateTCellStateMarkers:
    def test_restricts_to_available_genes(self):
        available = {"CCR7", "TCF7", "NOT_A_REAL_GENE"}
        result = validate_t_cell_state_markers(available, {"Naive_CM": ["CCR7", "TCF7", "SELL"]})
        assert result["Naive_CM"] == ["CCR7", "TCF7"]

    def test_raises_if_too_few_markers_present(self):
        available = {"CCR7"}
        with pytest.raises(PipelineError, match="only 1 marker"):
            validate_t_cell_state_markers(available, {"Naive_CM": ["CCR7", "TCF7", "SELL"]})


def _make_labeled_adata(n_per_state=100, n_genes=10, rng_seed=0):
    """Both states get their own elevated marker genes (state A: 3-5, state
    B: 0-2) -- a flat/no-signature "background" class would break Pearson-
    correlation classification (correlating against a near-zero-variance
    centroid is inherently unstable), which is not representative of the
    real use case: every curated T-cell state has its own distinguishing
    markers by construction."""
    rng = np.random.default_rng(rng_seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    state_a = rng.poisson(2, size=(n_per_state, n_genes)).astype(np.float32)
    state_b = rng.poisson(2, size=(n_per_state, n_genes)).astype(np.float32)
    state_a[:, 3:6] += 10  # genes 3-5 elevated in state A
    state_b[:, :3] += 10  # genes 0-2 elevated in state B
    X = np.vstack([state_a, state_b])
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell{i}" for i in range(X.shape[0])]
    labels = pd.Series(["A"] * n_per_state + ["B"] * n_per_state, index=adata.obs_names)
    return adata, labels


class TestComputeStateCentroids:
    def test_produces_one_row_per_label(self):
        adata, labels = _make_labeled_adata()
        centroids = compute_state_centroids(adata, labels)
        assert set(centroids.index) == {"A", "B"}
        assert centroids.shape[1] == adata.n_vars

    def test_respects_gene_subset(self):
        adata, labels = _make_labeled_adata()
        centroids = compute_state_centroids(adata, labels, gene_subset=["GENE0", "GENE1"])
        assert list(centroids.columns) == ["GENE0", "GENE1"]


class TestClassifyByNearestCentroid:
    def test_recovers_known_labels(self):
        adata, labels = _make_labeled_adata(rng_seed=1)
        centroids = compute_state_centroids(adata, labels)
        query = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
        result = classify_by_nearest_centroid(query, centroids)
        accuracy = (result["predicted_state"].to_numpy() == labels.to_numpy()).mean()
        assert accuracy > 0.9

    def test_recovers_labels_across_a_platform_scale_mismatch(self):
        """Regression test for a real bug found on real data: classifying a
        query dataset at a systematically different depth/scale than the
        reference produced a degenerate result (65% of all Xenium cells
        assigned to one state) when raw log1p-normalised values were
        compared directly, because log1p is not linear -- a fixed count
        -scale mismatch (e.g. reference target_sum=1e4 vs query
        target_sum=197) becomes a shape-distorting difference after
        log-transformation, not one Pearson correlation's row-centering
        alone can remove. The per-gene standardisation in
        classify_by_nearest_centroid must correct for this. Both reference
        and query are log1p'd here (as the real pipeline does), at
        deliberately different depths, to faithfully reproduce the failure
        mode -- not just a raw linear count-scale difference, which
        wouldn't exercise the actual bug."""
        rng = np.random.default_rng(3)
        n_per_state, n_genes = 150, 10
        ref_a = rng.poisson(50, size=(n_per_state, n_genes)).astype(np.float32)
        ref_b = rng.poisson(50, size=(n_per_state, n_genes)).astype(np.float32)
        ref_a[:, 3:6] += 200  # reference "depth": counts in the hundreds, like scRNA
        ref_b[:, :3] += 200
        ref_X = np.log1p(np.vstack([ref_a, ref_b]))
        ref_adata = ad.AnnData(X=ref_X)
        ref_adata.var_names = [f"GENE{i}" for i in range(n_genes)]
        ref_adata.obs_names = [f"cell{i}" for i in range(ref_X.shape[0])]
        ref_labels = pd.Series(["A"] * n_per_state + ["B"] * n_per_state, index=ref_adata.obs_names)
        centroids = compute_state_centroids(ref_adata, ref_labels)

        # Query "platform": same relative biological pattern, but ~25x
        # lower depth (single-digit counts, like a targeted spatial panel)
        # -- log1p'd separately, at its own low scale.
        rng2 = np.random.default_rng(4)
        query_a = rng2.poisson(2, size=(n_per_state, n_genes)).astype(np.float32)
        query_b = rng2.poisson(2, size=(n_per_state, n_genes)).astype(np.float32)
        query_a[:, 3:6] += 8
        query_b[:, :3] += 8
        query_X = np.log1p(np.vstack([query_a, query_b]))
        query = pd.DataFrame(
            query_X,
            columns=ref_adata.var_names,
            index=[f"lowdepth_cell{i}" for i in range(n_per_state * 2)],
        )
        query_labels = pd.Series(["A"] * n_per_state + ["B"] * n_per_state, index=query.index)

        result = classify_by_nearest_centroid(query, centroids)
        accuracy = (result["predicted_state"].to_numpy() == query_labels.to_numpy()).mean()
        assert accuracy > 0.8

    def test_raises_on_no_shared_genes(self):
        adata, labels = _make_labeled_adata()
        centroids = compute_state_centroids(adata, labels)
        query = pd.DataFrame(
            np.zeros((2, 2)), index=["c1", "c2"], columns=["UNRELATED1", "UNRELATED2"]
        )
        with pytest.raises(PipelineError, match="No shared genes"):
            classify_by_nearest_centroid(query, centroids)


class TestBenchmarkTransferConfidenceDegradation:
    def test_panel_with_discriminating_genes_classifies_well(self):
        """NOT a test that panel-restriction always hurts accuracy relative
        to a broader gene set: that would encode the same false intuition
        that caused a real bug on real data -- a nearest-
        centroid Pearson-correlation classifier weighs all included genes
        equally, so adding many uninformative genes can and does *dilute*
        signal and hurt accuracy relative to a well-chosen smaller gene set,
        confirmed on real data (the Xenium-panel-restricted set outperformed
        the ~38,600-gene whole reference before this was corrected to use a
        bounded HVG-based "broad" set instead). This test instead checks
        that a panel actually containing the true discriminating genes
        classifies well, and reports both accuracies without asserting a
        false ordering between them."""
        adata, labels = _make_labeled_adata(n_per_state=200, rng_seed=2)
        panel_genes = [
            "GENE0",
            "GENE1",
            "GENE2",
            "GENE3",
            "GENE4",
            "GENE5",
        ]  # the true discriminating genes
        result = benchmark_transfer_confidence_degradation(adata, labels, panel_genes, rng_seed=2)
        assert result["panel_restricted_accuracy"] > 0.8
        assert result["n_panel_overlap_genes"] == 6
        assert 0.0 <= result["broad_gene_set_accuracy"] <= 1.0


class TestLoadReferenceMatrix:
    def test_loads_a_synthetic_10x_style_matrix(self, tmp_path):
        n_cells, n_genes = 5, 4
        X = sparse.random(n_genes, n_cells, density=0.5, format="coo", random_state=0)
        X.data = np.round(X.data * 10).astype(np.float32) + 1

        with gzip.open(tmp_path / "matrix.mtx.gz", "wb") as f:
            mmwrite(f, X)

        features = pd.DataFrame(
            {
                "gene_id": [f"ENSG{i}" for i in range(n_genes)] + ["ENSG_AB"],
                "gene_symbol": [f"GENE{i}" for i in range(n_genes)] + ["CD3_hashtag"],
                "feature_type": ["Gene Expression"] * n_genes + ["Antibody Capture"],
            }
        )
        # matrix only has n_genes rows, so drop the antibody row to keep dims consistent for this test
        features = features.iloc[:n_genes].copy()
        features.to_csv(
            tmp_path / "features.tsv.gz", sep="\t", header=False, index=False, compression="gzip"
        )

        barcodes = pd.DataFrame({"barcode": [f"BC{i}-1" for i in range(n_cells)]})
        barcodes.to_csv(
            tmp_path / "barcodes.tsv.gz", sep="\t", header=False, index=False, compression="gzip"
        )

        adata = load_reference_matrix(tmp_path)
        assert adata.shape == (n_cells, n_genes)
        assert list(adata.var_names) == [f"GENE{i}" for i in range(n_genes)]

    def test_raises_on_missing_files(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            load_reference_matrix(tmp_path)
