"""Unit tests for xenium_tcr_ecology.external_checkpoint.program_transfer (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.external_checkpoint.program_transfer import (
    compute_mean_pairwise_correlation,
    compute_module_coherence,
    filter_expressed_genes,
)


class TestFilterExpressedGenes:
    def test_removes_zero_variance_genes(self):
        expr = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [0.0, 0.0, 0.0], "C": [5.0, 1.0, 9.0]})
        result = filter_expressed_genes(expr, ["A", "B", "C"])
        assert result == ["A", "C"]

    def test_keeps_all_genes_when_all_expressed(self):
        expr = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [4.0, 1.0, 7.0]})
        result = filter_expressed_genes(expr, ["A", "B"])
        assert set(result) == {"A", "B"}


class TestComputeMeanPairwiseCorrelation:
    def test_perfectly_correlated_genes_give_correlation_one(self):
        base = np.arange(50, dtype=float)
        expr = pd.DataFrame({"A": base, "B": base * 2, "C": base * 3})
        result = compute_mean_pairwise_correlation(expr, ["A", "B", "C"])
        assert result == pytest.approx(1.0)

    def test_uncorrelated_random_genes_give_near_zero(self):
        rng = np.random.default_rng(0)
        expr = pd.DataFrame({f"gene_{i}": rng.normal(size=200) for i in range(5)})
        result = compute_mean_pairwise_correlation(expr, list(expr.columns))
        assert abs(result) < 0.15

    def test_fewer_than_two_genes_gives_nan(self):
        expr = pd.DataFrame({"A": [1.0, 2.0, 3.0]})
        result = compute_mean_pairwise_correlation(expr, ["A"])
        assert np.isnan(result)


class TestComputeModuleCoherence:
    def test_real_coherent_module_beats_random_null(self):
        rng_data = np.random.default_rng(1)
        n = 300
        latent = rng_data.normal(size=n)
        # 4 genes strongly driven by a shared latent factor -- a real coherent module.
        module_genes = {
            f"module_{i}": latent + rng_data.normal(scale=0.2, size=n) for i in range(4)
        }
        # 30 independent noise genes for the null pool.
        noise_genes = {f"noise_{i}": rng_data.normal(size=n) for i in range(30)}
        expr = pd.DataFrame({**module_genes, **noise_genes})

        result = compute_module_coherence(
            expr,
            list(module_genes.keys()),
            list(noise_genes.keys()),
            np.random.default_rng(2),
            n_permutations=100,
        )
        assert result["observed_coherence"] > result["null_mean"]
        assert result["pvalue"] < 0.05

    def test_random_genes_do_not_beat_their_own_null(self):
        rng_data = np.random.default_rng(3)
        n = 300
        expr = pd.DataFrame({f"gene_{i}": rng_data.normal(size=n) for i in range(20)})
        test_genes = [f"gene_{i}" for i in range(4)]
        pool_genes = [f"gene_{i}" for i in range(4, 20)]

        result = compute_module_coherence(
            expr, test_genes, pool_genes, np.random.default_rng(4), n_permutations=100
        )
        assert result["pvalue"] > 0.05
