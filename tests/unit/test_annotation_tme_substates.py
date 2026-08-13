"""Unit tests for xenium_tcr_ecology.annotation.tme_substates (`06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.tme_substates import (
    SUBSTATE_MARKERS,
    compute_substate_scores,
    validate_substate_markers,
)


class TestValidateSubstateMarkers:
    def test_all_curated_substates_have_at_least_one_marker(self):
        for lineage, substates in SUBSTATE_MARKERS.items():
            for substate, genes in substates.items():
                assert len(genes) >= 1, f"{lineage}/{substate}"

    def test_restricts_to_available_genes(self):
        # Every curated gene present except one dropped from Macrophage's
        # (6-gene) set -- enough remain that no substate hits the 0-markers
        # error, so this isolates the restriction logic itself.
        available = {
            g for lineage in SUBSTATE_MARKERS.values() for genes in lineage.values() for g in genes
        }
        available.discard("TREM2")
        result = validate_substate_markers(available)
        assert "TREM2" not in result["Myeloid"]["Macrophage"]
        assert set(result["Myeloid"]["Macrophage"]) == set(
            SUBSTATE_MARKERS["Myeloid"]["Macrophage"]
        ) - {"TREM2"}

    def test_raises_if_a_substate_loses_all_markers(self):
        available = {"UNRELATED_GENE"}
        with pytest.raises(PipelineError, match="0 markers"):
            validate_substate_markers(available)


def _make_adata(n_cells=100, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    all_genes = sorted(
        {
            g
            for substates in SUBSTATE_MARKERS.values()
            for genes in substates.values()
            for g in genes
        }
    )
    genes = all_genes + [f"FILLER{i}" for i in range(60)]
    X = rng.poisson(3, size=(n_cells, len(genes))).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.layers["lognorm"] = np.log1p(X)
    return adata


class TestComputeSubstateScores:
    def test_returns_one_column_per_substate(self):
        adata = _make_adata()
        scores = compute_substate_scores(adata, layer="lognorm", gene_pool=list(adata.var_names))
        expected_cols = {
            f"{lineage}__{substate}_score"
            for lineage, substates in SUBSTATE_MARKERS.items()
            for substate in substates
        }
        assert set(scores.columns) == expected_cols
        assert len(scores) == adata.n_obs

    def test_cell_with_elevated_macrophage_markers_scores_higher_on_macrophage(self):
        adata = _make_adata(rng_seed=1)
        mac_genes = SUBSTATE_MARKERS["Myeloid"]["Macrophage"]
        mac_idx = [adata.var_names.get_loc(g) for g in mac_genes]

        X = adata.layers["lognorm"].copy()
        X[1, mac_idx] = X[1, mac_idx] + 5.0
        adata.layers["lognorm"] = X

        scores = compute_substate_scores(adata, layer="lognorm", gene_pool=list(adata.var_names))
        assert (
            scores.iloc[1]["Myeloid__Macrophage_score"]
            > scores.iloc[0]["Myeloid__Macrophage_score"]
        )
