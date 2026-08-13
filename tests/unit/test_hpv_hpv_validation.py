"""Unit tests for xenium_tcr_ecology.hpv.hpv_validation (`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.hpv.hpv_validation import (
    PROBE_POSITIVE_FRACTION_THRESHOLD,
    classify_validated_hpv_status,
    compute_fraction_probe_positive,
)


class TestComputeFractionProbePositive:
    def test_real_fraction_counts_any_gene_positive_per_cell(self):
        # 4 cells x 2 genes: cell 0 has E6>0, cell 1 has E7>0, cells 2-3 have neither.
        counts = np.array([[3, 0], [0, 2], [0, 0], [0, 0]])
        assert compute_fraction_probe_positive(counts) == 0.5

    def test_real_all_zero_gives_zero_fraction(self):
        counts = np.zeros((10, 2))
        assert compute_fraction_probe_positive(counts) == 0.0

    def test_real_all_positive_gives_one_fraction(self):
        counts = np.ones((5, 2))
        assert compute_fraction_probe_positive(counts) == 1.0


class TestClassifyValidatedHpvStatus:
    def test_real_confirmed_positive(self):
        result = classify_validated_hpv_status("Positive", True, 0.7)
        assert result == "confirmed_positive"

    def test_real_discordant_clinical_positive_probe_negative(self):
        result = classify_validated_hpv_status("Positive", True, 0.002)
        assert result == "discordant_clinical_positive_probe_negative"

    def test_real_clinical_positive_no_coverage(self):
        result = classify_validated_hpv_status("Positive", False, None)
        assert result == "clinical_positive_no_molecular_verification"

    def test_real_confirmed_negative(self):
        result = classify_validated_hpv_status("Negative", True, 0.0)
        assert result == "confirmed_negative"

    def test_real_discordant_clinical_negative_probe_positive(self):
        result = classify_validated_hpv_status("Negative", True, 0.3)
        assert result == "discordant_clinical_negative_probe_positive"

    def test_real_confirmed_negative_no_molecular_verification(self):
        result = classify_validated_hpv_status("Negative", False, None)
        assert result == "confirmed_negative_no_molecular_verification"

    def test_real_probe_positive_clinically_untested(self):
        result = classify_validated_hpv_status("Not Tested", True, 0.19)
        assert result == "probe_positive_clinically_untested"

    def test_real_probe_negative_clinically_untested(self):
        result = classify_validated_hpv_status("Not Tested", True, 0.0)
        assert result == "probe_negative_clinically_untested"

    def test_real_presumed_negative_unverifiable(self):
        result = classify_validated_hpv_status("Not Tested", False, None)
        assert result == "presumed_negative_unverifiable"

    def test_real_boundary_at_threshold_counts_as_positive(self):
        result = classify_validated_hpv_status("Positive", True, PROBE_POSITIVE_FRACTION_THRESHOLD)
        assert result == "confirmed_positive"
