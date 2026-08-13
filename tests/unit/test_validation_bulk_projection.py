"""Unit tests for xenium_tcr_ecology.validation.bulk_projection (`16_external_validation_and_generalisation/04_validate_ecosystem_signatures_in_bulk.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.validation.bulk_projection import (
    ECOSYSTEM_SIGNATURE_GENES,
    compute_signature_immune_correlation,
)


class TestEcosystemSignatureGenes:
    def test_real_mixed_non_specific_niche_is_excluded(self):
        assert "Mixed/non-specific niche" not in ECOSYSTEM_SIGNATURE_GENES

    def test_real_five_testable_ecosystems_present(self):
        assert len(ECOSYSTEM_SIGNATURE_GENES) == 5

    def test_real_immune_proxy_gene_not_in_any_signature(self):
        for genes in ECOSYSTEM_SIGNATURE_GENES.values():
            assert "PTPRC" not in genes


class TestComputeSignatureImmuneCorrelation:
    def test_real_positive_correlation_is_detected(self):
        rng = np.random.default_rng(0)
        n = 100
        immune_proxy = rng.normal(0, 1, n)
        signature_score = immune_proxy * 2 + rng.normal(0, 0.1, n)
        result = compute_signature_immune_correlation(
            pd.Series(signature_score), pd.Series(immune_proxy)
        )
        assert result["rho"] > 0.8
        assert result["pvalue"] < 0.05

    def test_real_no_correlation_gives_non_significant_result(self):
        rng = np.random.default_rng(1)
        n = 100
        signature_score = pd.Series(rng.normal(0, 1, n))
        immune_proxy = pd.Series(rng.normal(0, 1, n))
        result = compute_signature_immune_correlation(signature_score, immune_proxy)
        assert result["pvalue"] > 0.05
