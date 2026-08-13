"""Unit tests for xenium_tcr_ecology.validation.evidence_matrix (`16_external_validation_and_generalisation/07_generate_evidence_matrix.py`)."""

from __future__ import annotations

import pytest

from xenium_tcr_ecology.validation.evidence_matrix import (
    CLAIM_EVIDENCE_MATRIX,
    VALID_GRADES,
    validate_evidence_rows,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError


class TestValidateEvidenceRows:
    def test_real_well_formed_rows_pass(self):
        rows = [
            {
                "claim_id": "x",
                "discovery_evidence": "some",
                "sensitivity_evidence": "some",
                "replicate_evidence": "some",
                "external_validation_evidence": "some",
                "overall_evidence_grade": "strong",
            }
        ]
        validate_evidence_rows(rows)  # should not raise

    def test_real_invalid_grade_raises(self):
        rows = [
            {
                "claim_id": "x",
                "discovery_evidence": "some",
                "sensitivity_evidence": "some",
                "replicate_evidence": "some",
                "external_validation_evidence": "some",
                "overall_evidence_grade": "not_a_real_grade",
            }
        ]
        with pytest.raises(PipelineError):
            validate_evidence_rows(rows)

    def test_real_missing_field_raises(self):
        rows = [
            {
                "claim_id": "x",
                "discovery_evidence": "",
                "sensitivity_evidence": "some",
                "replicate_evidence": "some",
                "external_validation_evidence": "some",
                "overall_evidence_grade": "strong",
            }
        ]
        with pytest.raises(PipelineError):
            validate_evidence_rows(rows)


class TestClaimEvidenceMatrixConstant:
    def test_real_module_level_rows_are_well_formed(self):
        validate_evidence_rows(CLAIM_EVIDENCE_MATRIX)  # should not raise

    def test_real_every_registered_analysis_is_covered(self):
        import csv
        from pathlib import Path

        registry_path = Path(__file__).resolve().parents[2] / "governance" / "analysis_registry.tsv"
        with open(registry_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            registered_ids = {row["analysis_id"] for row in reader}
        covered_ids = {row["claim_id"] for row in CLAIM_EVIDENCE_MATRIX}
        assert registered_ids == covered_ids

    def test_real_grades_are_all_valid(self):
        for row in CLAIM_EVIDENCE_MATRIX:
            assert row["overall_evidence_grade"] in VALID_GRADES
