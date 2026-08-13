"""Unit tests for xenium_tcr_ecology.release.public_data_release (`17_statistical_closure_and_release/06_build_public_data_release.py`)."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.release.public_data_release import (
    ALL_SOURCE_FILES,
    build_public_data_release,
)


def _touch_source_files(project_root):
    for rel_path, _subdir in ALL_SOURCE_FILES:
        path = project_root / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        if rel_path.endswith(".yaml"):
            path.write_text("primary_contrasts: []\n")
        elif rel_path.endswith(".json"):
            path.write_text(json.dumps({"primary_analysis_ids": ["q1_framework_generalisation"]}))
        elif rel_path.endswith(".sha256"):
            path.write_text("deadbeef  placeholder\n")
        elif rel_path.endswith(".xlsx"):
            path.write_bytes(b"placeholder-xlsx-bytes")
        else:
            pd.DataFrame({"a": [1, 2]}).to_csv(path, sep="\t", index=False)


class TestBuildPublicDataRelease:
    def test_real_first_run_creates_a_real_manifest_and_license(self, tmp_path):
        _touch_source_files(tmp_path)
        summary = build_public_data_release(tmp_path)
        release_dir = tmp_path / "release" / "data"
        assert (release_dir / "MANIFEST.json").exists()
        assert (release_dir / "LICENSE").exists()
        assert (release_dir / "README.md").exists()
        assert (release_dir / "checksums.sha256").exists()
        assert summary["n_files"] == len(ALL_SOURCE_FILES)
        assert summary["license"] == "CC-BY-4.0"

    def test_real_files_land_in_their_declared_subdirectories(self, tmp_path):
        _touch_source_files(tmp_path)
        build_public_data_release(tmp_path)
        release_dir = tmp_path / "release" / "data"
        assert (release_dir / "primary_results" / "hpv_primary_contrasts.yaml").exists()
        assert (release_dir / "tables" / "Table_1_sample_manifest.tsv").exists()
        assert (release_dir / "metadata" / "data_dictionary.xlsx").exists()

    def test_real_rerun_with_unchanged_inputs_succeeds(self, tmp_path):
        _touch_source_files(tmp_path)
        build_public_data_release(tmp_path)
        summary = build_public_data_release(tmp_path)  # real, idempotent re-run
        assert summary["n_files"] == len(ALL_SOURCE_FILES)

    def test_real_rerun_with_a_changed_upstream_file_raises(self, tmp_path):
        _touch_source_files(tmp_path)
        build_public_data_release(tmp_path)
        changed_path = tmp_path / "tables" / "Table_1_sample_manifest.tsv"
        pd.DataFrame({"a": [999]}).to_csv(changed_path, sep="\t", index=False)
        with pytest.raises(PipelineError):
            build_public_data_release(tmp_path)

    def test_real_missing_upstream_file_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_public_data_release(tmp_path)

    def test_real_excludes_note_names_cell_level_data(self, tmp_path):
        _touch_source_files(tmp_path)
        build_public_data_release(tmp_path)
        manifest = json.loads((tmp_path / "release" / "data" / "MANIFEST.json").read_text())
        assert "cell-level" in manifest["excludes"]
