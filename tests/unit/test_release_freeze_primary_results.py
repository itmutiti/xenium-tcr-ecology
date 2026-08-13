"""Unit tests for xenium_tcr_ecology.release.freeze_primary_results (`17_statistical_closure_and_release/00_freeze_primary_results.py`)."""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.release.freeze_primary_results import (
    PRIMARY_RESULT_FILES,
    build_primary_results_freeze,
)


def _touch_primary_files(project_root):
    for rel_path in PRIMARY_RESULT_FILES:
        path = project_root / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        if rel_path.endswith(".yaml"):
            path.write_text("primary_contrasts: []\n")
        else:
            pd.DataFrame({"a": [1, 2]}).to_parquet(path)


class TestBuildPrimaryResultsFreeze:
    def test_real_first_run_creates_a_real_manifest(self, tmp_path):
        _touch_primary_files(tmp_path)
        summary = build_primary_results_freeze(tmp_path)
        manifest_path = tmp_path / "data" / "releases" / "final_primary" / "MANIFEST.json"
        assert manifest_path.exists()
        assert summary["n_files"] == len(PRIMARY_RESULT_FILES)

    def test_real_rerun_with_unchanged_inputs_succeeds(self, tmp_path):
        _touch_primary_files(tmp_path)
        build_primary_results_freeze(tmp_path)
        summary = build_primary_results_freeze(tmp_path)  # real, idempotent re-run
        assert summary["n_files"] == len(PRIMARY_RESULT_FILES)

    def test_real_rerun_with_a_changed_upstream_file_raises(self, tmp_path):
        _touch_primary_files(tmp_path)
        build_primary_results_freeze(tmp_path)
        changed_path = tmp_path / PRIMARY_RESULT_FILES[0]
        pd.DataFrame({"a": [999]}).to_parquet(changed_path)
        with pytest.raises(PipelineError):
            build_primary_results_freeze(tmp_path)

    def test_real_missing_upstream_file_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_primary_results_freeze(tmp_path)
