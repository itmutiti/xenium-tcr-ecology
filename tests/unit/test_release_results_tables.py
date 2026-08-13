"""Unit tests for xenium_tcr_ecology.release.results_tables (`17_statistical_closure_and_release/04_generate_results_tables.py`)."""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.release.results_tables import TABLE_MANIFEST, build_results_tables


def _touch_source_tables(project_root):
    for entry in TABLE_MANIFEST:
        path = project_root / entry["source_path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame({"a": [1, 2]})
        if entry["source_path"].endswith(".parquet"):
            df.to_parquet(path)
        else:
            df.to_csv(path, sep="\t", index=False)


class TestBuildResultsTables:
    def test_real_all_tables_written_as_real_tsv(self, tmp_path):
        _touch_source_tables(tmp_path)
        summary = build_results_tables(tmp_path)
        assert summary["n_tables"] == len(TABLE_MANIFEST)
        tables_dir = tmp_path / "tables"
        assert (tables_dir / "Table_1_sample_manifest.tsv").exists()

    def test_real_missing_source_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_results_tables(tmp_path)
