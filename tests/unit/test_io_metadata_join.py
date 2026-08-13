"""Unit tests for xenium_tcr_ecology.io.metadata_join (`03_spatialdata_import/04_attach_clinical_and_technical_metadata.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.io.metadata_join import METADATA_FIELDS, attach_metadata


def _write_manifest(tmp_path, rows):
    path = tmp_path / "sample_manifest.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _write_anndata(tmp_path, section_id, n_cells=5):
    adata = ad.AnnData(X=np.ones((n_cells, 3)))
    adata.obs["section_id"] = section_id
    path = tmp_path / f"{section_id}.h5ad"
    adata.write_h5ad(path)
    return path


class TestAttachMetadata:
    def _manifest_row(self, section_id="P01_run1"):
        row = {"section_id": section_id, "patient_id": "P01"}
        for f in METADATA_FIELDS:
            if f not in row:
                row[f] = "x"
        return row

    def test_joins_metadata_successfully(self, tmp_path):
        h5ad_path = _write_anndata(tmp_path, "P01_run1")
        manifest_path = _write_manifest(tmp_path, [self._manifest_row()])

        summary = attach_metadata(h5ad_path, manifest_path, tmp_path / "out.h5ad")
        assert summary["patient_id"] == "P01"

        result = ad.read_h5ad(tmp_path / "out.h5ad")
        assert (result.obs["patient_id"] == "P01").all()
        for f in METADATA_FIELDS:
            assert f in result.obs.columns

    def test_boolean_fields_are_real_booleans_not_strings(self, tmp_path):
        """Regression test: csv.DictReader reads every TSV field as a raw
        string; without an explicit cast, "True"/"False" text was being
        assigned directly into obs, silently producing a category column of
        strings rather than a real bool column -- present in every
        downstream artifact until `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` needed to use
        is_technical_replicate as an actual Python boolean mask."""
        h5ad_path = _write_anndata(tmp_path, "P01_run1")
        row = self._manifest_row()
        row["is_technical_replicate"] = "True"
        row["included_in_primary_hnscc_cohort"] = "False"
        manifest_path = _write_manifest(tmp_path, [row])

        attach_metadata(h5ad_path, manifest_path, tmp_path / "out.h5ad")
        result = ad.read_h5ad(tmp_path / "out.h5ad")

        assert result.obs["is_technical_replicate"].dtype == bool
        assert bool(result.obs["is_technical_replicate"].iloc[0]) is True
        assert result.obs["included_in_primary_hnscc_cohort"].dtype == bool
        assert bool(result.obs["included_in_primary_hnscc_cohort"].iloc[0]) is False

    def test_raises_on_duplicate_section_id_in_manifest(self, tmp_path):
        h5ad_path = _write_anndata(tmp_path, "P01_run1")
        manifest_path = _write_manifest(tmp_path, [self._manifest_row(), self._manifest_row()])

        with pytest.raises(PipelineError, match="duplicate section_id"):
            attach_metadata(h5ad_path, manifest_path, tmp_path / "out.h5ad")

    def test_raises_on_section_not_in_manifest(self, tmp_path):
        h5ad_path = _write_anndata(tmp_path, "P99_run1")
        manifest_path = _write_manifest(tmp_path, [self._manifest_row(section_id="P01_run1")])

        with pytest.raises(PipelineError, match="not found in"):
            attach_metadata(h5ad_path, manifest_path, tmp_path / "out.h5ad")

    def test_raises_on_mixed_section_ids_in_one_file(self, tmp_path):
        adata = ad.AnnData(X=np.ones((4, 3)))
        adata.obs["section_id"] = ["P01_run1", "P01_run1", "P02_run1", "P02_run1"]
        h5ad_path = tmp_path / "mixed.h5ad"
        adata.write_h5ad(h5ad_path)
        manifest_path = _write_manifest(
            tmp_path, [self._manifest_row("P01_run1"), self._manifest_row("P02_run1")]
        )

        with pytest.raises(PipelineError, match="more than one section_id"):
            attach_metadata(h5ad_path, manifest_path, tmp_path / "out.h5ad")

    def test_raises_on_missing_anndata_file(self, tmp_path):
        manifest_path = _write_manifest(tmp_path, [self._manifest_row()])
        with pytest.raises(PipelineError, match="not found"):
            attach_metadata(tmp_path / "nope.h5ad", manifest_path, tmp_path / "out.h5ad")
