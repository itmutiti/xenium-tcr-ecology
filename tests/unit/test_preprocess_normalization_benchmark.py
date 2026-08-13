"""Unit tests for xenium_tcr_ecology.preprocess.normalization_benchmark (`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.normalization_benchmark import (
    METHODS,
    PRIMARY_NORMALIZATION_LAYER,
    apply_primary_normalization_layer_decision,
    build_normalization_benchmark_summary,
    compute_replicate_stability,
    compute_technical_noise_correlation,
)


def _make_analysis_ready_fixture(tmp_path, n_cells_per_section=20, n_genes=6, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    sections = ["P01_run1", "P01_run2", "S1_run1"]  # P01 replicated, S1 not
    patients = {"P01_run1": "P01", "P01_run2": "P01", "S1_run1": "S1"}
    is_replicate = {"P01_run1": True, "P01_run2": True, "S1_run1": False}

    obs_names, section_ids, patient_ids, replicate_flags = [], [], [], []
    X_rows = []
    for section in sections:
        base = rng.poisson(5, size=(n_cells_per_section, n_genes)).astype(np.float32)
        for i in range(n_cells_per_section):
            obs_names.append(f"{section}_cell{i}")
            section_ids.append(section)
            patient_ids.append(patients[section])
            replicate_flags.append(is_replicate[section])
        X_rows.append(base)
    X = np.vstack(X_rows)

    adata = ad.AnnData(X=X.copy())
    adata.obs_names = obs_names
    adata.obs["section_id"] = pd.Categorical(section_ids)
    adata.obs["patient_id"] = pd.Categorical(patient_ids)
    adata.obs["is_technical_replicate"] = replicate_flags
    adata.var_names = [f"GENE{i}" for i in range(n_genes)]
    adata.var["is_exposure_gene"] = [True] * (n_genes - 1) + [False]

    adata.layers["counts"] = X.copy()
    adata.layers["lognorm"] = np.log1p(X)
    adata.layers["pearson_residuals"] = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-6)
    adata.layers["detected"] = (X > 0).astype(np.float32)

    analysis_ready_path = tmp_path / "analysis_ready.h5ad"
    adata.write_h5ad(analysis_ready_path)

    exclusion_log = pd.DataFrame({"cell_id": obs_names, "section_id": section_ids, "qc_pass": True})
    exclusion_log_path = tmp_path / "exclusion_log.tsv"
    exclusion_log.to_csv(exclusion_log_path, sep="\t", index=False)

    cell_qc = pd.DataFrame(
        {"control_probe_ratio": rng.uniform(0, 0.05, size=len(obs_names))}, index=obs_names
    )
    cell_qc_metrics_path = tmp_path / "cell_qc_metrics.parquet"
    cell_qc.to_parquet(cell_qc_metrics_path)

    return analysis_ready_path, exclusion_log_path, cell_qc_metrics_path


class TestComputeReplicateStability:
    def test_returns_one_row_per_replicate_patient(self, tmp_path):
        analysis_ready_path, exclusion_log_path, _ = _make_analysis_ready_fixture(tmp_path)
        result = compute_replicate_stability(analysis_ready_path, exclusion_log_path)
        assert len(result) == 1
        assert result.iloc[0]["patient_id"] == "P01"

    def test_computes_a_column_per_method(self, tmp_path):
        analysis_ready_path, exclusion_log_path, _ = _make_analysis_ready_fixture(tmp_path)
        result = compute_replicate_stability(analysis_ready_path, exclusion_log_path)
        for method in METHODS:
            assert f"{method}_replicate_r" in result.columns
            assert -1.0 <= result.iloc[0][f"{method}_replicate_r"] <= 1.0

    def test_raises_on_missing_analysis_ready_file(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            compute_replicate_stability(tmp_path / "nope.h5ad", tmp_path / "exclusion_log.tsv")


class TestComputeTechnicalNoiseCorrelation:
    def test_returns_one_row_per_method(self, tmp_path):
        analysis_ready_path, _, cell_qc_metrics_path = _make_analysis_ready_fixture(tmp_path)
        result = compute_technical_noise_correlation(analysis_ready_path, cell_qc_metrics_path)
        assert sorted(result["method"]) == sorted(METHODS)
        assert result["spearman_rho_vs_control_probe_ratio"].between(-1.0, 1.0).all()


class TestApplyPrimaryNormalizationLayerDecision:
    def test_writes_decision_onto_uns(self, tmp_path):
        analysis_ready_path, _, _ = _make_analysis_ready_fixture(tmp_path)
        assert "primary_normalization_layer" not in ad.read_h5ad(analysis_ready_path).uns

        wrote = apply_primary_normalization_layer_decision(analysis_ready_path)

        assert wrote is True
        reloaded = ad.read_h5ad(analysis_ready_path)
        assert reloaded.uns["primary_normalization_layer"] == PRIMARY_NORMALIZATION_LAYER

    def test_is_idempotent(self, tmp_path):
        analysis_ready_path, _, _ = _make_analysis_ready_fixture(tmp_path)
        assert apply_primary_normalization_layer_decision(analysis_ready_path) is True
        assert apply_primary_normalization_layer_decision(analysis_ready_path) is False


class TestBuildNormalizationBenchmarkSummary:
    def test_writes_outputs_and_returns_summary(self, tmp_path):
        (tmp_path / "data" / "objects").mkdir(parents=True)
        (tmp_path / "data" / "derived").mkdir(parents=True)
        analysis_ready_path, exclusion_log_path, cell_qc_metrics_path = (
            _make_analysis_ready_fixture(tmp_path)
        )
        analysis_ready_path.rename(tmp_path / "data" / "objects" / "analysis_ready.h5ad")
        exclusion_log_path.rename(tmp_path / "data" / "derived" / "exclusion_log.tsv")
        cell_qc_metrics_path.rename(tmp_path / "data" / "derived" / "cell_qc_metrics.parquet")

        output_path = tmp_path / "reports" / "preprocess" / "normalisation_benchmark"
        summary = build_normalization_benchmark_summary(tmp_path, output_path)

        assert summary["n_replicate_pairs"] == 1
        assert set(summary["median_replicate_r_by_method"].keys()) == set(METHODS)
        assert set(summary["abs_technical_noise_rho_by_method"].keys()) == set(METHODS)
        assert output_path.with_name(output_path.stem + "_replicate_stability.parquet").is_file()
        assert output_path.with_name(output_path.stem + "_technical_noise.parquet").is_file()
        reloaded = ad.read_h5ad(tmp_path / "data" / "objects" / "analysis_ready.h5ad")
        assert reloaded.uns["primary_normalization_layer"] == PRIMARY_NORMALIZATION_LAYER
