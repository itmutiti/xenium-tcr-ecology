"""Unit tests for xenium_tcr_ecology.governance (Project Setup and Governance charter, registry,
workflow init) and xenium_tcr_ecology.metadata (data dictionary)."""

from __future__ import annotations

import pytest
import yaml

from xenium_tcr_ecology.governance.charter import build_project_charter
from xenium_tcr_ecology.governance.registry import compile_analysis_registry
from xenium_tcr_ecology.governance.workflow_init import initialise_workflow
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.metadata.data_dictionary import compile_data_dictionary


class TestProjectCharter:
    def test_always_succeeds_and_records_primary_questions(self, tmp_path):
        (tmp_path / "governance").mkdir()
        output_path = tmp_path / "charter.yaml"
        charter = build_project_charter(tmp_path, output_path)
        assert set(charter["primary_questions"]) == {"Q1", "Q2", "Q3"}
        assert output_path.is_file()


class TestAnalysisRegistry:
    def _input(self, tmp_path, hpv_primary_count=1):
        analyses = [
            {
                "analysis_id": "a1",
                "phase": "10_niche_and_ecosystem_discovery",
                "hypothesis": "h1",
                "unit_of_analysis": "clone",
                "primary_endpoint": "e1",
                "multiplicity_family": "primary (Q2)",
            }
        ]
        for i in range(hpv_primary_count):
            analyses.append(
                {
                    "analysis_id": f"hpv_{i}",
                    "phase": "15_hpv_stratified_analysis",
                    "hypothesis": "h",
                    "unit_of_analysis": "patient",
                    "primary_endpoint": "e",
                    "multiplicity_family": "primary (HPV) contrast",
                }
            )
        path = tmp_path / "registry_input.yaml"
        path.write_text(yaml.dump({"analyses": analyses}))
        return path

    def test_compiles_within_hpv_cap(self, tmp_path):
        input_path = self._input(tmp_path, hpv_primary_count=2)
        summary = compile_analysis_registry(
            input_path,
            tmp_path / "out.tsv",
            project_root=tmp_path,
            registered_by="tester",
            registered_date="2026-01-01",
        )
        assert summary["analyses_registered"] == 3
        assert summary["hpv_primary_contrasts_reserved"] == 2

    def test_raises_over_hpv_cap(self, tmp_path):
        input_path = self._input(tmp_path, hpv_primary_count=3)
        with pytest.raises(PipelineError, match="exceeding the cap"):
            compile_analysis_registry(
                input_path,
                tmp_path / "out.tsv",
                project_root=tmp_path,
                registered_by="tester",
                registered_date="2026-01-01",
            )


class TestWorkflowInit:
    def test_raises_if_snakefile_missing(self, tmp_path):
        (tmp_path / "config").mkdir()
        (tmp_path / "config" / "config.yaml").write_text("x: 1")
        (tmp_path / "config" / "global_seed.yaml").write_text("default_seed: 42")
        with pytest.raises(PipelineError, match="Snakefile"):
            initialise_workflow(tmp_path)

    def test_writes_profile_once_then_skips(self, tmp_path):
        (tmp_path / "config").mkdir()
        (tmp_path / "Snakefile").write_text("rule all:\n    input: []\n")
        (tmp_path / "config" / "config.yaml").write_text("x: 1")
        (tmp_path / "config" / "global_seed.yaml").write_text("default_seed: 42")

        summary1 = initialise_workflow(tmp_path)
        assert summary1["profile_written"] is True
        assert summary1["default_seed"] == 42

        summary2 = initialise_workflow(tmp_path)
        assert summary2["profile_written"] is False  # already present, not forced


class TestDataDictionary:
    def test_compiles_and_writes_xlsx(self, tmp_path):
        pytest.importorskip("openpyxl")
        data = {
            "tables": [
                {
                    "table_name": "test_table.tsv",
                    "fields": [
                        {
                            "field": "col_a",
                            "unit": "categorical",
                            "allowed_values": "x|y",
                            "missingness_code": "none",
                            "derivation_rule": "raw input",
                        }
                    ],
                }
            ]
        }
        input_path = tmp_path / "dict_input.yaml"
        input_path.write_text(yaml.dump(data))
        output_path = tmp_path / "dict.xlsx"

        summary = compile_data_dictionary(input_path, output_path, project_root=tmp_path)
        assert summary["tables_documented"] == 1
        assert output_path.is_file()

    def test_raises_on_column_drift_from_real_table(self, tmp_path):
        """The whole point of the cross-check: if the real table on disk
        doesn't match what's documented, this must fail loudly, not
        silently document a stale schema."""
        pytest.importorskip("openpyxl")
        (tmp_path / "metadata").mkdir()
        (tmp_path / "metadata" / "drifted.tsv").write_text("col_a\tcol_new\nval\tval\n")

        data = {
            "tables": [
                {
                    "table_name": "drifted.tsv",
                    "fields": [
                        {
                            "field": "col_a",
                            "unit": "categorical",
                            "allowed_values": "x",
                            "missingness_code": "none",
                            "derivation_rule": "raw input",
                        }
                    ],
                }
            ]
        }
        input_path = tmp_path / "dict_input.yaml"
        input_path.write_text(yaml.dump(data))

        with pytest.raises(PipelineError, match="out of sync"):
            compile_data_dictionary(input_path, tmp_path / "dict.xlsx", project_root=tmp_path)


class TestRealProjectPhase1Outputs:
    def test_real_project_charter_records_primary_questions(self):
        """Exercises the real project_charter.yaml already written for this
        project."""
        project_root = find_project_root()
        charter_path = project_root / "project_charter.yaml"
        if not charter_path.is_file():
            pytest.skip("project_charter.yaml not yet generated")
        charter = yaml.safe_load(charter_path.read_text())
        assert set(charter["primary_questions"]) == {"Q1", "Q2", "Q3"}
