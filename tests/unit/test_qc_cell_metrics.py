"""Unit tests for xenium_tcr_ecology.qc.cell_metrics (`04_quality_control/00_compute_cell_level_qc_metrics.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.qc.cell_metrics import build_cell_qc_report, compute_cell_qc_metrics


def _make_adata(n_cells=6, n_genes=4):
    rng = np.random.default_rng(0)
    X = csr_matrix(rng.poisson(2, size=(n_cells, n_genes)).astype(float))
    adata = ad.AnnData(X=X)
    adata.obs["section_id"] = ["P01_run1"] * 3 + ["P02_run1"] * 3
    adata.obs["patient_id"] = ["P01"] * 3 + ["P02"] * 3
    # transcript_counts is gene-expression-only in real 10x cells.parquet data
    # (confirmed against real data, exactly equal to sum(X) for all 1,186,916
    # cells in this project) -- control_probe_counts/control_codeword_counts
    # are separate, sibling fields, not included in transcript_counts.
    adata.obs["transcript_counts"] = np.asarray(X.sum(axis=1)).ravel()
    adata.obs["control_probe_counts"] = 5
    adata.obs["control_codeword_counts"] = 5
    adata.obs["cell_area"] = 100.0
    adata.obs["nucleus_area"] = 40.0
    return adata


class TestComputeCellQcMetrics:
    def test_computes_expected_columns(self):
        adata = _make_adata()
        metrics = compute_cell_qc_metrics(adata)

        for col in [
            "n_genes_detected",
            "counts_from_expression_matrix",
            "counts_discrepancy",
            "control_probe_ratio",
            "control_codeword_ratio",
            "transcript_density_per_um2",
            "nucleus_to_cell_area_ratio",
        ]:
            assert col in metrics.columns

        # transcript_counts is gene-expression-only, so it must equal
        # counts_from_expression_matrix even when control_probe_counts /
        # control_codeword_counts are nonzero -- this is the real behaviour
        # confirmed against the full combined object (0 discrepancy across
        # all 1,186,916 cells), not an assumption.
        assert (metrics["counts_discrepancy"] == 0).all()
        assert np.allclose(metrics["nucleus_to_cell_area_ratio"], 0.4)

    def test_flags_a_genuine_export_integrity_error(self):
        """counts_discrepancy must still have teeth: if transcript_counts and
        sum(X) diverge (e.g. a cell/gene lost or duplicated during
        AnnData export), that must be visible, not silently absorbed."""
        adata = _make_adata()
        adata.obs["transcript_counts"] = adata.obs["transcript_counts"] + 7
        metrics = compute_cell_qc_metrics(adata)
        assert (metrics["counts_discrepancy"] == 7).all()

    def test_raises_on_missing_required_column(self):
        adata = ad.AnnData(X=np.ones((3, 2)))
        with pytest.raises(PipelineError, match="missing required column"):
            compute_cell_qc_metrics(adata)


class TestBuildCellQcReport:
    def test_writes_report_and_summary(self, tmp_path):
        adata = _make_adata()
        h5ad_path = tmp_path / "combined.h5ad"
        adata.write_h5ad(h5ad_path)

        summary = build_cell_qc_report(h5ad_path, tmp_path / "cell_qc_metrics.parquet")
        assert summary["n_cells"] == 6
        assert summary["n_sections"] == 2
        assert (tmp_path / "cell_qc_metrics.parquet").is_file()

    def test_raises_on_missing_input(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            build_cell_qc_report(tmp_path / "nope.h5ad", tmp_path / "out.parquet")
