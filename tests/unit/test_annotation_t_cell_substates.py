"""Unit tests for xenium_tcr_ecology.annotation.t_cell_substates (`06_cell_type_annotation/04_resolve_t_cell_substates.R`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.t_cell_substates import compute_treg_score


def _make_adata(n_cells=100, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    genes = ["FOXP3", "IL2RA", "CTLA4", "CD4", "CD8A"] + [f"FILLER{i}" for i in range(45)]
    X = rng.poisson(3, size=(n_cells, len(genes))).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.layers["lognorm"] = np.log1p(X)
    return adata


class TestComputeTregScore:
    def test_returns_a_score_per_cell(self):
        adata = _make_adata()
        result = compute_treg_score(adata, layer="lognorm", gene_pool=list(adata.var_names))
        assert len(result) == adata.n_obs
        assert np.isfinite(result.to_numpy()).all()

    def test_cell_with_elevated_treg_markers_scores_higher(self):
        adata = _make_adata(rng_seed=1)
        treg_idx = [adata.var_names.get_loc(g) for g in ["FOXP3", "IL2RA", "CTLA4"]]
        X = adata.layers["lognorm"].copy()
        X[1, treg_idx] = X[1, treg_idx] + 5.0
        adata.layers["lognorm"] = X

        result = compute_treg_score(adata, layer="lognorm", gene_pool=list(adata.var_names))
        assert result.iloc[1] > result.iloc[0]

    def test_raises_if_too_few_treg_markers_present(self):
        adata = ad.AnnData(X=np.ones((5, 2), dtype=np.float32))
        adata.var_names = ["FOXP3", "OTHER"]
        adata.obs_names = [f"c{i}" for i in range(5)]
        adata.layers["lognorm"] = adata.X.copy()
        with pytest.raises(PipelineError, match="Treg marker"):
            compute_treg_score(adata, layer="lognorm", gene_pool=list(adata.var_names))
