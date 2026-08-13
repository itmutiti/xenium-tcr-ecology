"""Unit tests for xenium_tcr_ecology.hpv.claim_strength (`15_hpv_stratified_analysis/06_prepare_hpv_claim_strength_table.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.hpv.claim_strength import grade_group_comparison


class TestGradeGroupComparison:
    def test_real_stress_tested_stable_direction_and_bootstrap_excludes_zero_is_exploratory(self):
        stress_info = {
            "bootstrap_ci_low": 0.01,
            "bootstrap_ci_high": 0.1,
            "n_lopo_direction_flips": 0,
        }
        assert grade_group_comparison(stress_info) == "exploratory"

    def test_real_stress_tested_bootstrap_crosses_zero_is_unsuitable(self):
        stress_info = {
            "bootstrap_ci_low": -0.5,
            "bootstrap_ci_high": 0.3,
            "n_lopo_direction_flips": 0,
        }
        assert grade_group_comparison(stress_info) == "unsuitable_for_inference"

    def test_real_stress_tested_lopo_flips_is_unsuitable_even_if_bootstrap_excludes_zero(self):
        stress_info = {
            "bootstrap_ci_low": 0.01,
            "bootstrap_ci_high": 0.1,
            "n_lopo_direction_flips": 2,
        }
        assert grade_group_comparison(stress_info) == "unsuitable_for_inference"

    def test_real_never_stress_tested_is_unsuitable(self):
        assert grade_group_comparison(None) == "unsuitable_for_inference"

    def test_real_group_comparison_is_never_graded_supported(self):
        stress_info = {
            "bootstrap_ci_low": 0.5,
            "bootstrap_ci_high": 0.9,
            "n_lopo_direction_flips": 0,
        }
        assert grade_group_comparison(stress_info) != "supported"
