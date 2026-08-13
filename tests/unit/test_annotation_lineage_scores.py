"""Unit tests for xenium_tcr_ecology.annotation.lineage_scores (`06_cell_type_annotation/02_score_major_lineages.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.lineage_scores import (
    compute_lineage_scores,
    load_major_lineage_gene_sets,
)


class TestLoadMajorLineageGeneSets:
    def test_parses_semicolon_joined_markers(self):
        registry = pd.DataFrame(
            {
                "cell_identity": ["T_cell", "B_cell", "Lymphatic_endothelial"],
                "hierarchy_level": ["major_lineage", "major_lineage", "substate"],
                "markers": ["CD3D;CD3E", "CD19;MS4A1", "PROX1;LYVE1"],
            }
        )
        result = load_major_lineage_gene_sets(registry)
        assert result == {"T_cell": ["CD3D", "CD3E"], "B_cell": ["CD19", "MS4A1"]}
        assert "Lymphatic_endothelial" not in result  # substate excluded

    def test_raises_if_no_major_lineage_rows(self):
        registry = pd.DataFrame(
            {"cell_identity": ["X"], "hierarchy_level": ["substate"], "markers": ["A;B"]}
        )
        with pytest.raises(PipelineError, match="no major_lineage"):
            load_major_lineage_gene_sets(registry)

    def test_raises_if_a_lineage_has_too_few_markers(self):
        registry = pd.DataFrame(
            {"cell_identity": ["X"], "hierarchy_level": ["major_lineage"], "markers": ["A"]}
        )
        with pytest.raises(PipelineError, match="below the minimum"):
            load_major_lineage_gene_sets(registry)


def _make_adata(n_cells=100, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    genes = ["CD3D", "CD3E", "CD2", "CD19", "MS4A1", "CD79A"] + [f"FILLER{i}" for i in range(30)]
    X = rng.poisson(3, size=(n_cells, len(genes))).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.layers["lognorm"] = np.log1p(X)
    return adata


class TestComputeLineageScores:
    def test_returns_one_column_per_lineage(self):
        adata = _make_adata()
        gene_sets = {"T_cell": ["CD3D", "CD3E", "CD2"], "B_cell": ["CD19", "MS4A1", "CD79A"]}
        scores = compute_lineage_scores(
            adata, layer="lognorm", gene_sets=gene_sets, gene_pool=list(adata.var_names)
        )
        assert set(scores.columns) == {"T_cell_lineage_score", "B_cell_lineage_score"}
        assert len(scores) == adata.n_obs

    def test_cell_expressing_only_t_cell_genes_scores_higher_on_t_cell(self):
        adata = _make_adata(n_cells=50, rng_seed=1)
        gene_sets = {"T_cell": ["CD3D", "CD3E", "CD2"], "B_cell": ["CD19", "MS4A1", "CD79A"]}
        t_idx = [adata.var_names.get_loc(g) for g in gene_sets["T_cell"]]

        X = adata.layers["lognorm"].copy()
        X[1, t_idx] = X[1, t_idx] + 5.0
        adata.layers["lognorm"] = X

        scores = compute_lineage_scores(
            adata, layer="lognorm", gene_sets=gene_sets, gene_pool=list(adata.var_names)
        )
        assert scores.iloc[1]["T_cell_lineage_score"] > scores.iloc[0]["T_cell_lineage_score"]

    def test_raises_if_lineage_genes_missing_from_var_names(self):
        adata = _make_adata()
        gene_sets = {"Nonsense": ["NOT_A_GENE"]}
        with pytest.raises(PipelineError, match="below the minimum"):
            compute_lineage_scores(
                adata, layer="lognorm", gene_sets=gene_sets, gene_pool=list(adata.var_names)
            )
