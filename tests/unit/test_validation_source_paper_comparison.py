"""Unit tests for xenium_tcr_ecology.validation.source_paper_comparison (`16_external_validation_and_generalisation/06_compare_with_source_paper_results.py`)."""

from __future__ import annotations

from pathlib import Path

import pytest

from xenium_tcr_ecology.validation.source_paper_comparison import (
    SOURCE_PAPER_COMPARISON,
    VALID_STATUSES,
    build_source_paper_comparison,
    validate_comparison_rows,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError


class TestValidateComparisonRows:
    def test_real_well_formed_rows_pass(self):
        rows = [
            {
                "figure_or_table": "Figure X",
                "comparison_status": "reproduced",
                "this_project_reference": "Project Setup and Governance",
                "note": "some note",
                "novel_contribution": "something new",
            },
        ]
        validate_comparison_rows(rows)  # should not raise

    def test_real_invalid_status_raises(self):
        rows = [
            {
                "figure_or_table": "Figure X",
                "comparison_status": "not_a_real_status",
                "this_project_reference": "Project Setup and Governance",
                "note": "some note",
                "novel_contribution": "n/a",
            },
        ]
        with pytest.raises(PipelineError):
            validate_comparison_rows(rows)

    def test_real_missing_field_raises(self):
        rows = [
            {
                "figure_or_table": "Figure X",
                "comparison_status": "reproduced",
                "this_project_reference": "",
                "note": "some note",
                "novel_contribution": "n/a",
            },
        ]
        with pytest.raises(PipelineError):
            validate_comparison_rows(rows)


class TestSourcePaperComparisonConstant:
    def test_real_module_level_rows_are_well_formed(self):
        validate_comparison_rows(SOURCE_PAPER_COMPARISON)  # should not raise

    def test_real_statuses_are_all_valid(self):
        for row in SOURCE_PAPER_COMPARISON:
            assert row["comparison_status"] in VALID_STATUSES


class TestBuildSourcePaperComparison:
    def test_succeeds_and_writes_full_table(self, tmp_path):
        summary = build_source_paper_comparison(tmp_path)

        assert summary["n_claims_compared"] == len(SOURCE_PAPER_COMPARISON)
        output_path = Path(summary["output_path"])
        assert output_path.is_file()
        import pandas as pd

        written = pd.read_csv(output_path, sep="\t")
        assert len(written) == len(SOURCE_PAPER_COMPARISON)
        assert set(written["figure_or_table"]) == {
            row["figure_or_table"] for row in SOURCE_PAPER_COMPARISON
        }
