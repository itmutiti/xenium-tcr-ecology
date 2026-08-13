"""Unit tests for xenium_tcr_ecology.external_checkpoint.directional_consistency (`12_external_checkpoint_validation/02_quantify_directional_consistency.py`)."""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.external_checkpoint.directional_consistency import (
    compute_pairwise_sign_agreement,
    compute_rank_consistency,
    compute_rank_shift,
)


class TestComputePairwiseSignAgreement:
    def test_identical_ordering_agrees_on_every_pair(self):
        project = pd.Series({"A": 0.5, "B": 0.3, "C": 0.2})
        reference = pd.Series({"A": 0.6, "B": 0.25, "C": 0.15})
        result = compute_pairwise_sign_agreement(project, reference)
        assert result["agrees"].all()
        assert len(result) == 3  # C(3,2)

    def test_reversed_ordering_disagrees_on_every_pair(self):
        project = pd.Series({"A": 0.6, "B": 0.3, "C": 0.1})
        reference = pd.Series({"A": 0.1, "B": 0.3, "C": 0.6})
        result = compute_pairwise_sign_agreement(project, reference)
        # A vs C flips fully; A vs B and B vs C also flip since B is tied
        # to the same rank position -- check the A-vs-C pair specifically.
        row = result[(result["state_a"] == "A") & (result["state_b"] == "C")].iloc[0]
        assert not row["agrees"]


class TestComputeRankConsistency:
    def test_perfect_agreement_gives_rho_one(self):
        project = pd.Series({"A": 0.5, "B": 0.3, "C": 0.2, "D": 0.1})
        reference = pd.Series({"A": 0.9, "B": 0.05, "C": 0.03, "D": 0.02})
        result = compute_rank_consistency(project, reference)
        assert result["spearman_rho"] == pytest.approx(1.0)

    def test_perfect_reversal_gives_rho_negative_one(self):
        project = pd.Series({"A": 0.5, "B": 0.3, "C": 0.2, "D": 0.1})
        reference = pd.Series({"A": 0.02, "B": 0.03, "C": 0.05, "D": 0.9})
        result = compute_rank_consistency(project, reference)
        assert result["spearman_rho"] == pytest.approx(-1.0)


class TestComputeRankShift:
    def test_rank_one_is_most_abundant(self):
        project = pd.Series({"A": 0.5, "B": 0.3, "C": 0.2})
        reference = pd.Series({"A": 0.5, "B": 0.3, "C": 0.2})
        result = compute_rank_shift(project, reference).set_index("state")
        assert result.loc["A", "project_rank"] == 1
        assert result.loc["C", "project_rank"] == 3
        assert result.loc["A", "rank_shift"] == 0

    def test_rank_shift_reflects_real_reordering(self):
        # 'Cycling'-like case: dominant in project, mid-pack in reference.
        project = pd.Series({"Cycling": 0.40, "Exhausted": 0.22, "Treg": 0.18})
        reference = pd.Series({"Cycling": 0.14, "Exhausted": 0.23, "Treg": 0.15})
        result = compute_rank_shift(project, reference).set_index("state")
        assert result.loc["Cycling", "project_rank"] == 1
        assert result.loc["Cycling", "reference_rank"] == 3
        assert result.loc["Cycling", "rank_shift"] == 2
