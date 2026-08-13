"""Unit tests for xenium_tcr_ecology.qc.apply_filters (`04_quality_control/07_apply_qc_filters_with_audit_trail.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.qc.apply_filters import (
    EXCLUSION_REASON_COLUMNS,
    _local_cell_id,
    apply_qc_filters,
    evaluate_thresholds,
)

PROFILE = {
    "min_transcript_counts": 10,
    "min_genes_detected": 3,
    "max_control_probe_ratio": 0.05,
    "max_control_codeword_ratio": 0.05,
    "section_relative_min_counts_z": -3.5,
}


def _make_cell_qc():
    return pd.DataFrame(
        {
            "transcript_counts": [50, 3, 50, 50, 50],
            "n_genes_detected": [10, 10, 1, 10, 10],
            "control_probe_ratio": [0.0, 0.0, 0.0, 0.20, 0.0],
            "control_codeword_ratio": [0.0, 0.0, 0.0, 0.0, 0.0],
            "z_counts": [0.0, 0.0, 0.0, 0.0, -10.0],
        },
        index=["clean", "low_count", "low_genes", "high_ctrl", "z_outlier"],
    )


class TestEvaluateThresholds:
    def test_flags_low_transcript_count(self):
        reasons = evaluate_thresholds(_make_cell_qc(), PROFILE)
        assert reasons.loc["low_count", "low_transcript_count"]

    def test_flags_low_genes_detected(self):
        reasons = evaluate_thresholds(_make_cell_qc(), PROFILE)
        assert reasons.loc["low_genes", "low_genes_detected"]

    def test_flags_high_control_probe_ratio(self):
        reasons = evaluate_thresholds(_make_cell_qc(), PROFILE)
        assert reasons.loc["high_ctrl", "high_control_probe_ratio"]

    def test_flags_section_relative_outlier(self):
        reasons = evaluate_thresholds(_make_cell_qc(), PROFILE)
        assert reasons.loc["z_outlier", "section_relative_low_count_outlier"]

    def test_does_not_flag_a_clean_cell(self):
        reasons = evaluate_thresholds(_make_cell_qc(), PROFILE)
        assert not reasons.loc["clean"].any()


class TestLocalCellId:
    def test_strips_section_prefix(self):
        index = pd.Index(["P01_run1_aaadggoi-1", "P09_run2_zzzzzzzz-1"])
        section_ids = pd.Series(["P01_run1", "P09_run2"], index=index)
        result = _local_cell_id(index, section_ids)
        assert result.tolist() == ["aaadggoi-1", "zzzzzzzz-1"]


def _make_combined_adata():
    rng = np.random.default_rng(0)
    n_cells = 6
    X = csr_matrix(rng.poisson(2, size=(n_cells, 3)).astype(float))
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.obs["patient_id"] = ["P01"] * 3 + ["P02"] * 3
    return adata


class TestApplyQcFilters:
    def test_removes_only_excluded_cells(self, tmp_path):
        adata = _make_combined_adata()
        h5ad_path = tmp_path / "combined.h5ad"
        adata.write_h5ad(h5ad_path)

        exclusion_log = pd.DataFrame(
            {"qc_pass": [True, True, False, True, True, True]}, index=adata.obs_names
        )
        out_path = tmp_path / "qc_filtered.h5ad"
        filtered = apply_qc_filters(h5ad_path, exclusion_log, out_path)

        assert filtered.n_obs == 5
        assert "cell2" not in filtered.obs_names
        assert out_path.is_file()

    def test_raises_if_an_entire_patient_would_be_excluded(self, tmp_path):
        adata = _make_combined_adata()
        h5ad_path = tmp_path / "combined.h5ad"
        adata.write_h5ad(h5ad_path)

        # P02 is cells 3, 4, 5 -- exclude all of them
        exclusion_log = pd.DataFrame(
            {"qc_pass": [True, True, True, False, False, False]}, index=adata.obs_names
        )
        with pytest.raises(PipelineError, match="P02"):
            apply_qc_filters(h5ad_path, exclusion_log, tmp_path / "qc_filtered.h5ad")

    def test_raises_if_exclusion_log_missing_cells(self, tmp_path):
        adata = _make_combined_adata()
        h5ad_path = tmp_path / "combined.h5ad"
        adata.write_h5ad(h5ad_path)

        exclusion_log = pd.DataFrame({"qc_pass": [True, True]}, index=["cell0", "cell1"])
        with pytest.raises(PipelineError, match="exclusion log"):
            apply_qc_filters(h5ad_path, exclusion_log, tmp_path / "qc_filtered.h5ad")


def test_exclusion_reason_columns_are_consistent_set():
    assert len(EXCLUSION_REASON_COLUMNS) == len(set(EXCLUSION_REASON_COLUMNS))
