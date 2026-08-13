"""Unit tests for xenium_tcr_ecology.preprocess.release_freeze (`05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`)."""

from __future__ import annotations

import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.release_freeze import (
    compute_file_hash,
    freeze_primary_analysis_matrix,
    merge_program_scores_into_obs,
)


class TestComputeFileHash:
    def test_matches_known_sha256(self, tmp_path):
        path = tmp_path / "f.txt"
        path.write_bytes(b"hello world")
        expected = hashlib.sha256(b"hello world").hexdigest()
        assert compute_file_hash(path) == expected

    def test_different_content_gives_different_hash(self, tmp_path):
        p1, p2 = tmp_path / "a.txt", tmp_path / "b.txt"
        p1.write_bytes(b"content a")
        p2.write_bytes(b"content b")
        assert compute_file_hash(p1) != compute_file_hash(p2)


class TestMergeProgramScoresIntoObs:
    def test_adds_score_columns(self):
        adata = ad.AnnData(X=np.ones((3, 2)))
        adata.obs_names = ["c1", "c2", "c3"]
        scores = pd.DataFrame({"cytotoxicity_score": [0.1, 0.2, 0.3]}, index=["c1", "c2", "c3"])
        result = merge_program_scores_into_obs(adata, scores)
        assert list(result.obs["cytotoxicity_score"]) == [0.1, 0.2, 0.3]

    def test_raises_if_a_cell_is_missing_from_scores(self):
        adata = ad.AnnData(X=np.ones((3, 2)))
        adata.obs_names = ["c1", "c2", "c3"]
        scores = pd.DataFrame({"cytotoxicity_score": [0.1, 0.2]}, index=["c1", "c2"])
        with pytest.raises(PipelineError, match="no entry"):
            merge_program_scores_into_obs(adata, scores)

    def test_raises_on_column_name_collision(self):
        adata = ad.AnnData(X=np.ones((2, 2)))
        adata.obs_names = ["c1", "c2"]
        adata.obs["cytotoxicity_score"] = [1.0, 2.0]
        scores = pd.DataFrame({"cytotoxicity_score": [0.1, 0.2]}, index=["c1", "c2"])
        with pytest.raises(PipelineError, match="collide"):
            merge_program_scores_into_obs(adata, scores)


class TestFreezePrimaryAnalysisMatrix:
    def _make_project(self, tmp_path):
        (tmp_path / "data" / "objects").mkdir(parents=True)
        (tmp_path / "data" / "derived").mkdir(parents=True)

        adata = ad.AnnData(X=np.ones((4, 2), dtype=np.float32))
        adata.obs_names = [f"c{i}" for i in range(4)]
        adata.var_names = ["GENE1", "GENE2"]
        adata.layers["lognorm"] = np.log1p(adata.X)
        adata.uns["primary_normalization_layer"] = "lognorm"
        adata.write_h5ad(tmp_path / "data" / "objects" / "analysis_ready.h5ad")

        scores = pd.DataFrame({"cytotoxicity_score": [0.1, 0.2, 0.3, 0.4]}, index=adata.obs_names)
        scores.to_parquet(tmp_path / "data" / "derived" / "program_scores.parquet")
        return tmp_path

    def test_writes_matrix_manifest_and_checksums(self, tmp_path):
        project_root = self._make_project(tmp_path)
        release_dir = tmp_path / "data" / "releases" / "v1_primary_analysis"

        summary = freeze_primary_analysis_matrix(project_root, release_dir)

        assert summary["n_cells"] == 4
        assert (release_dir / "primary_analysis_matrix.h5ad").is_file()
        assert (release_dir / "MANIFEST.json").is_file()
        assert (release_dir / "checksums.sha256").is_file()

        manifest = json.loads((release_dir / "MANIFEST.json").read_text())
        assert manifest["primary_normalization_layer"] == "lognorm"
        assert "cytotoxicity_score" in manifest["obs_columns"]

    def test_manifest_hash_matches_actual_file(self, tmp_path):
        project_root = self._make_project(tmp_path)
        release_dir = tmp_path / "data" / "releases" / "v1_primary_analysis"
        summary = freeze_primary_analysis_matrix(project_root, release_dir)

        actual_hash = compute_file_hash(release_dir / "primary_analysis_matrix.h5ad")
        assert summary["matrix_hash"] == actual_hash

    def test_raises_on_missing_analysis_ready(self, tmp_path):
        (tmp_path / "data" / "derived").mkdir(parents=True)
        with pytest.raises(PipelineError, match="not found"):
            freeze_primary_analysis_matrix(tmp_path, tmp_path / "release")
