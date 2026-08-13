"""Unit tests for xenium_tcr_ecology.validation.plan (`16_external_validation_and_generalisation/00_define_validation_claims.py`)."""

from __future__ import annotations

import pytest

from xenium_tcr_ecology.validation.plan import VALIDATION_CLAIMS, validate_claims_well_formed
from xenium_tcr_ecology.infra.exceptions import PipelineError


class TestValidateClaimsWellFormed:
    def test_real_well_formed_claims_pass(self):
        claims = [
            {
                "claim_id": "a",
                "claim": "some claim",
                "validation_dataset": "some dataset",
                "validation_method": "some method",
                "success_criterion": "some criterion",
                "phase_reference": "16.01",
            },
        ]
        validate_claims_well_formed(claims)  # should not raise

    def test_real_missing_field_raises(self):
        claims = [
            {
                "claim_id": "a",
                "claim": "some claim",
                "validation_dataset": "some dataset",
                "validation_method": "some method",
                "success_criterion": "",
                "phase_reference": "16.01",
            },
        ]
        with pytest.raises(PipelineError):
            validate_claims_well_formed(claims)

    def test_real_duplicate_claim_id_raises(self):
        claims = [
            {
                "claim_id": "a",
                "claim": "x",
                "validation_dataset": "y",
                "validation_method": "z",
                "success_criterion": "w",
                "phase_reference": "16.01",
            },
            {
                "claim_id": "a",
                "claim": "x2",
                "validation_dataset": "y2",
                "validation_method": "z2",
                "success_criterion": "w2",
                "phase_reference": "16.02",
            },
        ]
        with pytest.raises(PipelineError):
            validate_claims_well_formed(claims)


class TestRealClaimsConstant:
    def test_real_module_level_claims_are_well_formed(self):
        validate_claims_well_formed(VALIDATION_CLAIMS)  # should not raise

    def test_real_q1_framework_generalisation_claim_is_present(self):
        assert any(c["claim_id"] == "q1_framework_generalisation" for c in VALIDATION_CLAIMS)
