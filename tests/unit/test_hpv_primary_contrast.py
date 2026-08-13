"""Unit tests for xenium_tcr_ecology.hpv.primary_contrast (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`)."""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.hpv.primary_contrast import (
    assign_contrast_group,
    build_hpv_contrast_config,
    build_primary_hpv_contrasts,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError


class TestAssignContrastGroup:
    @pytest.mark.parametrize(
        "status,expected",
        [
            ("confirmed_positive", "positive"),
            ("probe_positive_clinically_untested", "positive"),
            ("confirmed_negative", "negative"),
            ("confirmed_negative_no_molecular_verification", "negative"),
            ("presumed_negative_unverifiable", "negative"),
            ("probe_negative_clinically_untested", "negative"),
            ("discordant_clinical_positive_probe_negative", "excluded"),
            ("discordant_clinical_negative_probe_positive", "excluded"),
            ("clinical_positive_no_molecular_verification", "excluded"),
        ],
    )
    def test_real_every_phase_15_00_status_maps_to_a_group(self, status, expected):
        assert assign_contrast_group(status) == expected


class TestBuildHpvContrastConfig:
    def test_real_patients_are_grouped_correctly(self):
        hpv_status = pd.DataFrame(
            {
                "patient_id": ["P09", "P01", "P20"],
                "validated_hpv_status": [
                    "confirmed_positive",
                    "discordant_clinical_positive_probe_negative",
                    "confirmed_negative",
                ],
            }
        )
        result = build_hpv_contrast_config(hpv_status, registered_date="2026-07-11")
        assert result["positive_group"]["patient_ids"] == ["P09"]
        assert result["negative_group"]["patient_ids"] == ["P20"]
        assert result["excluded_patients"]["patient_ids"] == ["P01"]
        assert result["n_positive"] == 1
        assert result["n_negative"] == 1
        assert result["n_excluded"] == 1
        assert result["registered_date"] == "2026-07-11"
        assert result["model"] is None
        assert result["minimum_detectable_effect"] is None

    def test_real_empty_positive_group_raises(self):
        hpv_status = pd.DataFrame(
            {
                "patient_id": ["P01"],
                "validated_hpv_status": ["discordant_clinical_positive_probe_negative"],
            }
        )
        with pytest.raises(PipelineError):
            build_hpv_contrast_config(hpv_status, registered_date="2026-07-11")

    def test_real_empty_negative_group_raises(self):
        hpv_status = pd.DataFrame(
            {
                "patient_id": ["P09"],
                "validated_hpv_status": ["confirmed_positive"],
            }
        )
        with pytest.raises(PipelineError):
            build_hpv_contrast_config(hpv_status, registered_date="2026-07-11")


class TestBuildPrimaryHpvContrasts:
    def test_real_cap_blocks_a_second_different_contrast(self, tmp_path):
        (tmp_path / "metadata").mkdir()
        hpv_status_path = tmp_path / "metadata" / "hpv_status_validated.tsv"
        pd.DataFrame(
            {
                "patient_id": ["P09", "P20"],
                "validated_hpv_status": ["confirmed_positive", "confirmed_negative"],
            }
        ).to_csv(hpv_status_path, sep="\t", index=False)

        (tmp_path / "governance").mkdir()
        existing_path = tmp_path / "governance" / "hpv_primary_contrasts.yaml"
        existing_path.write_text(
            "primary_contrasts:\n"
            "  - contrast_id: some_other_contrast\n"
            "    registered_date: '2026-01-01'\n"
            "    model: null\n"
            "    minimum_detectable_effect: null\n"
            "  - contrast_id: yet_another_contrast\n"
            "    registered_date: '2026-01-02'\n"
            "    model: null\n"
            "    minimum_detectable_effect: null\n"
        )

        with pytest.raises(PipelineError):
            build_primary_hpv_contrasts(tmp_path)

    def test_real_rerun_is_idempotent(self, tmp_path):
        (tmp_path / "metadata").mkdir()
        hpv_status_path = tmp_path / "metadata" / "hpv_status_validated.tsv"
        pd.DataFrame(
            {
                "patient_id": ["P09", "P20"],
                "validated_hpv_status": ["confirmed_positive", "confirmed_negative"],
            }
        ).to_csv(hpv_status_path, sep="\t", index=False)

        first = build_primary_hpv_contrasts(tmp_path)
        second = build_primary_hpv_contrasts(tmp_path)
        assert first == second
