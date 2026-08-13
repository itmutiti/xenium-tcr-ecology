"""Unit tests for xenium_tcr_ecology.io.combine_sections (`03_spatialdata_import/05_build_combined_analysis_object.py`) and
xenium_tcr_ecology.io.r_export (`03_spatialdata_import/06_export_r_interoperability_objects.py`), using synthetic fixtures."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.io import mmread

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.io.combine_sections import build_combined_object
from xenium_tcr_ecology.io.r_export import export_section_r_interop


def _write_section_anndata(path, section_id, patient_id, cell_ids, gene_names):
    n_cells, n_genes = len(cell_ids), len(gene_names)
    adata = ad.AnnData(X=np.arange(n_cells * n_genes).reshape(n_cells, n_genes).astype(float))
    adata.var_names = gene_names
    adata.obs["cell_id"] = cell_ids
    adata.obs["section_id"] = section_id
    adata.obs["patient_id"] = patient_id
    adata.obsm["spatial"] = np.column_stack([np.arange(n_cells), np.arange(n_cells) * 2])
    adata.write_h5ad(path)
    return adata


class TestBuildCombinedObject:
    def test_combines_sections_with_globally_unique_ids(self, tmp_path):
        anndata_root = tmp_path / "anndata"
        anndata_root.mkdir()
        genes = ["G1", "G2"]
        # Same raw cell_id "cell_0" appears in both sections -- must not collide after combining.
        _write_section_anndata(
            anndata_root / "P01_run1.h5ad", "P01_run1", "P01", ["cell_0", "cell_1"], genes
        )
        _write_section_anndata(
            anndata_root / "P02_run1.h5ad", "P02_run1", "P02", ["cell_0", "cell_1"], genes
        )

        summary = build_combined_object(
            anndata_root, tmp_path / "combined.h5ad", tmp_path / "panel_membership.parquet"
        )
        assert summary["n_sections"] == 2
        assert summary["n_cells_total"] == 4
        assert summary["n_patients"] == 2
        assert summary["n_genes_union"] == 2
        assert summary["n_genes_core"] == 2

        combined = ad.read_h5ad(tmp_path / "combined.h5ad")
        assert not combined.obs_names.duplicated().any()
        assert "P01_run1_cell_0" in combined.obs_names
        assert "P02_run1_cell_0" in combined.obs_names

    def test_handles_mismatched_gene_panels_via_outer_join(self, tmp_path):
        """Gene panels legitimately differ across sections (patient-specific/
        batch-specific CDR3 probes, confirmed against real project data) --
        this must succeed via an outer join, not raise, and must record
        which genes were actually in each section's panel so a 0-fill is
        never later confused with 'measured and absent'."""
        anndata_root = tmp_path / "anndata"
        anndata_root.mkdir()
        _write_section_anndata(
            anndata_root / "P01_run1.h5ad", "P01_run1", "P01", ["c0"], ["G1", "G2"]
        )
        _write_section_anndata(
            anndata_root / "P02_run1.h5ad", "P02_run1", "P02", ["c0"], ["G1", "G3"]
        )

        panel_path = tmp_path / "panel_membership.parquet"
        summary = build_combined_object(anndata_root, tmp_path / "combined.h5ad", panel_path)
        assert summary["n_genes_union"] == 3  # G1, G2, G3
        assert summary["n_genes_core"] == 1  # only G1 shared by both

        combined = ad.read_h5ad(tmp_path / "combined.h5ad")
        assert set(combined.var_names) == {"G1", "G2", "G3"}

        membership = pd.read_parquet(panel_path)
        assert membership.loc["G2", "P01_run1"] and not membership.loc["G2", "P02_run1"]
        assert membership.loc["G3", "P02_run1"] and not membership.loc["G3", "P01_run1"]
        assert membership.loc["G1", "P01_run1"] and membership.loc["G1", "P02_run1"]

    def test_raises_on_missing_required_column(self, tmp_path):
        anndata_root = tmp_path / "anndata"
        anndata_root.mkdir()
        adata = ad.AnnData(X=np.ones((2, 2)))
        adata.write_h5ad(anndata_root / "bad.h5ad")

        with pytest.raises(PipelineError, match="missing required column"):
            build_combined_object(
                anndata_root, tmp_path / "combined.h5ad", tmp_path / "panel.parquet"
            )

    def test_raises_on_empty_directory(self, tmp_path):
        anndata_root = tmp_path / "anndata"
        anndata_root.mkdir()
        with pytest.raises(PipelineError, match="No .h5ad files"):
            build_combined_object(
                anndata_root, tmp_path / "combined.h5ad", tmp_path / "panel.parquet"
            )


class TestExportSectionRInterop:
    def test_exports_matrix_market_and_metadata(self, tmp_path):
        h5ad_path = tmp_path / "P01_run1.h5ad"
        _write_section_anndata(h5ad_path, "P01_run1", "P01", ["c0", "c1", "c2"], ["G1", "G2"])
        output_dir = tmp_path / "r_export"

        summary = export_section_r_interop(h5ad_path, output_dir)
        assert summary["n_cells"] == 3
        assert summary["n_genes"] == 2

        import gzip

        with gzip.open(output_dir / "matrix.mtx.gz", "rb") as f:
            matrix = mmread(f)
        # Matrix Market export is genes x cells (transposed from AnnData's cells x genes).
        assert matrix.shape == (2, 3)

        with gzip.open(output_dir / "barcodes.tsv.gz", "rt") as f:
            barcodes = f.read().strip().splitlines()
        assert len(barcodes) == 3

        obs = pd.read_parquet(output_dir / "obs_metadata.parquet")
        assert "x_centroid" in obs.columns and "y_centroid" in obs.columns

    def test_raises_on_missing_input(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            export_section_r_interop(tmp_path / "nope.h5ad", tmp_path / "out")
