"""Unit tests for xenium_tcr_ecology.preprocess.program_scores (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.program_scores import (
    MIN_GENES_PER_PROGRAM,
    PROGRAM_GENE_SETS,
    compute_program_scores,
    validate_gene_sets,
)


class TestValidateGeneSets:
    def test_keeps_only_genes_present_in_panel(self):
        available = {
            "GZMA",
            "GZMB",
            "GZMK",
            "PRF1",
            "GNLY",
            "NKG7",
            "KLRD1",
            "KLRB1",
            "KLRC1",
            "FGFBP2",
        }
        result = validate_gene_sets(
            available, {"cytotoxicity": PROGRAM_GENE_SETS["cytotoxicity"] + ["NOT_A_REAL_GENE"]}
        )
        assert "NOT_A_REAL_GENE" not in result["cytotoxicity"]
        assert set(result["cytotoxicity"]) == available

    def test_raises_if_too_few_genes_survive(self):
        available = {"GZMA"}
        with pytest.raises(PipelineError, match="below the minimum"):
            validate_gene_sets(available, {"cytotoxicity": ["GZMA", "NOT_PRESENT"]})

    def test_all_curated_programs_have_enough_real_panel_genes(self):
        """Every program's curated set, as actually shipped, must clear the
        minimum-genes bar against itself (i.e. is internally consistent) --
        catches a program accidentally curated down to <2 genes."""
        for program, genes in PROGRAM_GENE_SETS.items():
            assert len(genes) >= MIN_GENES_PER_PROGRAM, program


def _make_adata(n_cells=200, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    genes = sorted(set().union(*PROGRAM_GENE_SETS.values())) + [f"FILLER{i}" for i in range(50)]
    X = rng.poisson(3, size=(n_cells, len(genes))).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.layers["lognorm"] = np.log1p(X)
    return adata


class TestComputeProgramScores:
    def test_returns_one_score_column_per_program(self):
        adata = _make_adata()
        scores = compute_program_scores(adata, layer="lognorm", gene_pool=list(adata.var_names))
        assert set(scores.columns) == {f"{p}_score" for p in PROGRAM_GENE_SETS}
        assert len(scores) == adata.n_obs

    def test_scores_are_finite(self):
        adata = _make_adata()
        scores = compute_program_scores(adata, layer="lognorm", gene_pool=list(adata.var_names))
        assert np.isfinite(scores.to_numpy()).all()

    def test_cell_expressing_only_cytotoxicity_genes_scores_higher_on_cytotoxicity(self):
        """Sanity check on real behaviour, not just shape: a cell with
        strongly elevated cytotoxicity-gene expression (and nothing else)
        should score higher on cytotoxicity than a flat-expression cell."""
        adata = _make_adata(n_cells=100, rng_seed=1)
        cyto_genes = PROGRAM_GENE_SETS["cytotoxicity"]
        cyto_idx = [adata.var_names.get_loc(g) for g in cyto_genes]

        X = adata.layers["lognorm"].copy()
        # Cell 0: baseline. Cell 1: elevated cytotoxicity genes only.
        X[1, cyto_idx] = X[1, cyto_idx] + 5.0
        adata.layers["lognorm"] = X

        scores = compute_program_scores(adata, layer="lognorm", gene_pool=list(adata.var_names))
        assert scores.iloc[1]["cytotoxicity_score"] > scores.iloc[0]["cytotoxicity_score"]
