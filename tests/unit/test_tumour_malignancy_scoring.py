"""Unit tests for xenium_tcr_ecology.tumour.malignancy_scoring (`07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`)."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.malignancy_scoring import (
    MIN_MARKERS_PRESENT,
    _zscore_within_patient,
    combine_malignancy_evidence,
    compute_hpv_score,
    compute_patient_clonality_score,
    compute_tumour_marker_score,
)


class TestZscoreWithinPatient:
    def test_all_nan_group_stays_nan_not_zero(self):
        # A patient with no HPV panel at all must stay NaN through
        # standardisation, not be coerced to 0.0 ("exactly average") --
        # that would silently pull "no data" into the combined mean as if
        # it were real evidence of average malignancy.
        series = pd.Series([np.nan, np.nan, 1.0, 3.0])
        patients = pd.Series(["A", "A", "B", "B"])
        result = _zscore_within_patient(series, patients)
        assert result[:2].isna().all()
        assert result[2:].notna().all()

    def test_large_scale_patient_does_not_affect_small_scale_patient(self):
        # pandas' default (sample, ddof=1) std of a 2-point group [0, 2] is
        # sqrt(2), giving z = +-1/sqrt(2), not +-1.0.
        series = pd.Series([0.0, 2.0, 0.0, 200.0])
        patients = pd.Series(["A", "A", "B", "B"])
        result = _zscore_within_patient(series, patients)
        assert result[0] == pytest.approx(-1 / np.sqrt(2))
        assert result[1] == pytest.approx(1 / np.sqrt(2))


def _make_adata(n_cells=200, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    genes = (
        ["EPCAM", "MET", "ERBB2", "EGFR", "KRT7"]
        + [f"HPV16_{x}" for x in ("E1", "E2", "E5", "E6", "E7")]
        + [f"FILLER{i}" for i in range(40)]
    )
    X = rng.poisson(3, size=(n_cells, len(genes))).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.layers["lognorm"] = np.log1p(X)
    adata.obs["section_id"] = ["S1"] * (n_cells // 2) + ["S2"] * (n_cells - n_cells // 2)
    return adata


class TestComputeTumourMarkerScore:
    def test_raises_when_too_few_markers_present(self):
        adata = _make_adata()
        adata = adata[:, ["FILLER0", "FILLER1"]].copy()
        adata.layers["lognorm"] = np.log1p(adata.X)
        with pytest.raises(PipelineError, match="below minimum"):
            compute_tumour_marker_score(adata, "lognorm", gene_pool=list(adata.var_names))

    def test_returns_one_score_per_cell(self):
        adata = _make_adata()
        result = compute_tumour_marker_score(adata, "lognorm", gene_pool=list(adata.var_names))
        assert len(result) == adata.n_obs
        assert result.index.equals(adata.obs_names)
        assert MIN_MARKERS_PRESENT >= 2


class TestComputeHpvScore:
    def test_null_for_sections_without_hpv_panel(self):
        adata = _make_adata()
        panel_membership = pd.DataFrame(
            {"S1": [True] * 5, "S2": [False] * 5},
            index=[f"HPV16_{x}" for x in ("E1", "E2", "E5", "E6", "E7")],
        )
        result = compute_hpv_score(adata, "lognorm", list(adata.var_names), panel_membership)
        s1_cells = adata.obs_names[adata.obs["section_id"] == "S1"]
        s2_cells = adata.obs_names[adata.obs["section_id"] == "S2"]
        assert result.loc[s1_cells].notna().all()
        assert result.loc[s2_cells].isna().all()

    def test_works_when_gene_pool_excludes_hpv_genes(self):
        # Regression test: the real gene_pool passed by
        # build_malignancy_score_report is restricted to `biological_gene`
        # features only (`05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`'s composition-bias convention), which
        # deliberately excludes the HPV probes themselves. scanpy's
        # score_genes bins control genes using `gene_pool` alone, then
        # looks up each `gene_list` gene's own bin in that binning -- if
        # gene_list is not a subset of gene_pool, this raises a KeyError.
        # The first real run hit exactly this.
        adata = _make_adata()
        hpv_genes = [g for g in adata.var_names if g.startswith("HPV16_")]
        gene_pool_excluding_hpv = [g for g in adata.var_names if g not in hpv_genes]
        panel_membership = pd.DataFrame(
            {"S1": [True] * len(hpv_genes), "S2": [True] * len(hpv_genes)}, index=hpv_genes
        )
        result = compute_hpv_score(adata, "lognorm", gene_pool_excluding_hpv, panel_membership)
        assert result.notna().all()

    def test_all_null_when_no_hpv_genes_in_panel(self):
        adata = _make_adata()
        adata = adata[:, [g for g in adata.var_names if not g.startswith("HPV16_")]].copy()
        adata.layers["lognorm"] = np.log1p(adata.X)
        panel_membership = pd.DataFrame({"S1": [], "S2": []})
        result = compute_hpv_score(adata, "lognorm", list(adata.var_names), panel_membership)
        assert result.isna().all()


class TestComputePatientClonalityScore:
    def test_single_patient_cluster_scores_one(self):
        clusters = pd.Series(["0", "0", "1", "1"])
        patients = pd.Series(["P1", "P1", "P2", "P2"])
        result = compute_patient_clonality_score(clusters, patients)
        assert (result == 1.0).all()

    def test_mixed_cluster_scores_below_one(self):
        clusters = pd.Series(["0", "0", "0", "0"])
        patients = pd.Series(["P1", "P1", "P2", "P3"])
        result = compute_patient_clonality_score(clusters, patients)
        assert (result == 0.5).all()

    def test_returns_float_dtype_for_categorical_cluster_labels(self):
        # Regression test: Leiden cluster labels (the real input, from
        # scanpy) are categorical dtype. An earlier version returned a
        # Categorical Series from `.map()` in this case, which raised
        # TypeError in the combine step's `.mean()` on real data.
        clusters = pd.Series(["0", "0", "1", "1"], dtype="category")
        patients = pd.Series(["P1", "P1", "P2", "P2"])
        result = compute_patient_clonality_score(clusters, patients)
        assert pd.api.types.is_float_dtype(result)


class TestCombineMalignancyEvidence:
    def test_hpv_null_cells_still_get_a_combined_score(self):
        idx = [f"c{i}" for i in range(10)]
        tumour = pd.Series(np.linspace(-1, 1, 10), index=idx)
        hpv = pd.Series([np.nan] * 5 + list(np.linspace(-1, 1, 5)), index=idx)
        emt = pd.Series(np.linspace(1, -1, 10), index=idx)
        clonality = pd.Series([0.5] * 10, index=idx)
        patients = pd.Series(["P1"] * 10, index=idx)
        result = combine_malignancy_evidence(tumour, hpv, emt, clonality, patients)
        assert result["malignancy_score"].notna().all()

    def test_malignancy_probability_is_percentile_rank_in_0_1(self):
        idx = [f"c{i}" for i in range(20)]
        rng = np.random.default_rng(1)
        tumour = pd.Series(rng.normal(size=20), index=idx)
        hpv = pd.Series(rng.normal(size=20), index=idx)
        emt = pd.Series(rng.normal(size=20), index=idx)
        clonality = pd.Series(rng.uniform(0.2, 1.0, size=20), index=idx)
        patients = pd.Series((["P1"] * 10) + (["P2"] * 10), index=idx)
        result = combine_malignancy_evidence(tumour, hpv, emt, clonality, patients)
        assert result["malignancy_probability"].min() > 0
        assert result["malignancy_probability"].max() <= 1.0

    def test_higher_raw_scores_get_higher_malignancy_probability(self):
        idx = ["low", "high"]
        tumour = pd.Series([-5.0, 5.0], index=idx)
        hpv = pd.Series([-5.0, 5.0], index=idx)
        emt = pd.Series([-5.0, 5.0], index=idx)
        clonality = pd.Series([0.2, 1.0], index=idx)
        patients = pd.Series(["P1", "P1"], index=idx)
        result = combine_malignancy_evidence(tumour, hpv, emt, clonality, patients)
        assert (
            result.loc["high", "malignancy_probability"]
            > result.loc["low", "malignancy_probability"]
        )

    def test_one_high_scale_patient_does_not_invert_another_patients_ranking(self):
        # Regression test: this is the exact real-data failure mode. Patient A's HPV signal is
        # 0-2 (a real, modest malignant-vs-normal contrast); patient B's is
        # 0-200 (a much larger absolute scale, unrelated to per-cell
        # malignancy). A pooled/global z-score lets patient B dominate the
        # mean/std and pushes patient A's higher-than-its-own-
        # baseline cell toward a LOW value; a within-patient z-score must
        # not let this happen.
        idx = ["a_low", "a_high", "b_low", "b_high"]
        hpv = pd.Series([0.0, 2.0, 0.0, 200.0], index=idx)
        tumour = pd.Series([0.0, 0.0, 0.0, 0.0], index=idx)
        emt = pd.Series([0.0, 0.0, 0.0, 0.0], index=idx)
        clonality = pd.Series([0.5, 0.5, 0.5, 0.5], index=idx)
        patients = pd.Series(["A", "A", "B", "B"], index=idx)
        result = combine_malignancy_evidence(tumour, hpv, emt, clonality, patients)
        # Patient A's own higher-HPV cell must score above patient A's own
        # lower-HPV cell -- a pooled z-score would instead make "a_high"
        # (raw value 2.0) look tiny next to "b_high" (raw value 200.0) and
        # could distort A's within-patient ordering via a shared mean/std.
        assert result.loc["a_high", "hpv_score"] > result.loc["a_low", "hpv_score"]
        assert result.loc["b_high", "hpv_score"] > result.loc["b_low", "hpv_score"]
        # And patient A's own within-patient standardised value should
        # depend only on A's own 2-point spread (+-1/sqrt(2)), not be
        # compressed toward 0 by patient B's much larger scale.
        assert result.loc["a_high", "hpv_score"] == pytest.approx(1 / np.sqrt(2))
        assert result.loc["a_low", "hpv_score"] == pytest.approx(-1 / np.sqrt(2))
