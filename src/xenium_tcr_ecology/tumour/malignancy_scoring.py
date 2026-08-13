"""Malignancy probability estimation (`07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`).

Combines four evidence sources into a within-epithelial-subset malignancy
score, all computed on `07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py`'s `epithelial_subset.h5ad` (i.e.
"malignant vs normal epithelium," a within-lineage question, not
"epithelial vs other lineage," which was already `06_cell_type_annotation/02_score_major_lineages.py`'s job):

1. **Tumour marker score** (`tumour_marker_score`): `scanpy.tl.score_genes`
   over `06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s `Epithelial_Tumour` registry markers (EPCAM, MET,
   ERBB2, EGFR, KRT7), against the same 399-gene `biological_gene` control
   pool used throughout this project (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`, 6.02). Documented
   limitation, carried over from `marker_registry.py`: none of the
   canonical squamous-differentiation markers (TP63, KRT5, KRT14, KRT17,
   KRT19, KRT8, KRT18, SFN) are in this panel, so this captures general
   epithelial/oncogene signal (ERBB2/EGFR/MET overexpression is an
   established HNSCC malignancy signal), not confirmed
   squamous-malignancy-specificity.
2. **HPV score** (`hpv_score`): `scanpy.tl.score_genes` over the 8
   HPV16_* probes, computed only for cells whose section's panel
   physically includes them (`results/tables/03_spatialdata_import/gene_panel_membership.parquet`
   -- confirmed that HPV probe presence is panel/patient-specific,
   like the CDR3 probes, not a universal zero-vs-nonzero read on every
   section). NaN, not 0, for cells from HPV-probe-free sections: absence
   of an HPV signal where the probe was never physically present is
   uninformative about malignancy, not evidence against it, and must not
   silently penalise those cells in the combined score.
3. **EMT/stress score** (`emt_stress_score`): mean of `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s
   pre-computed `emt_score` and `stress_score` program scores (already in
   `primary_analysis_matrix.h5ad.obs`, carried through `07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py`'s
   subsetting) -- both established in the HNSCC/carcinoma literature as
   elevated in malignant/invasive-front epithelium relative to normal
   epithelium.
4. **Patient-clonality score** (`patient_clonality_score`): substitutes the
   specification's "reference mapping" component. This project's
   only whole-transcriptome external reference (GSE287301, `06_cell_type_annotation/03_map_external_scrna_reference.py`) is
   T-cell-only (366,632 T cells, no epithelial/malignant cells at all --
   confirmed against that module's own docstring), so there is no
   available reference dataset to map epithelial cells against. Acquiring
   a new external malignant-epithelium reference (e.g. Puram et al. 2017)
   is a substantial new data-acquisition undertaking, out of scope for
   this milestone -- not silently approximated with partial/synthetic
   data here. In its place: each cell's own `07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py` joint cluster's
   patient-exclusivity fraction (the per-cluster dominant-patient
   fraction, `epithelial_subset.compute_joint_cluster_patient_dominance`,
   assigned per-cell rather than aggregated) is used as a principled,
   in-scope proxy for the same underlying biology reference-mapping-based
   malignancy scoring is usually built on: patient-specific copy-number/
   genotype drives both (a) failure to map cleanly onto a shared reference
   and (b) failure to cluster jointly with other patients' cells, so a
   cell's own cluster's patient-exclusivity is evidence of the same
   underlying process, not an unrelated stand-in.

All four raw scores are z-scored within each patient (not pooled across the
whole epithelial subset) before combination. This is stronger than the
"within the epithelial subset, not the whole dataset" scoping alone would
give, and was itself a bug found and fixed on the data: raw HPV16_* transcript counts vary
by 2-3 orders of magnitude across clinically HPV+ patients for reasons
unrelated to per-cell malignancy (viral copy number, tumour purity) --
confirmed (median 38 counts/epithelial-cell in one patient, 0 in
several other clinically HPV+ patients). A single pooled z-score let the
one high-signal patient dominate the mean/std, silently inverting the
expected direction for every other HPV+ patient. This is the same
principle already applied for the same reason in `06_cell_type_annotation/03_map_external_scrna_reference.py` (per-platform
standardisation before cross-platform comparison) and is applied to all
four components here, not only the one that exposed the bug, since the
underlying question ("which of this patient's epithelial cells look more
malignant than that patient's own normal-epithelium baseline") is
inherently within-patient. `malignancy_score` is the row-wise mean of
whichever within-patient-z-scored components are available for a cell
(only `hpv_score` is ever missing); `malignancy_probability` is its
percentile rank within the epithelial
subset, matching this project's established convention for turning a raw
combined score into an interpretable [0, 1] value (`06_cell_type_annotation/06_integrate_annotation_evidence.py`'s
`marker_margin`).

This phase deliberately outputs a continuous probability, not a hard
malignant/normal call: `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py` ("Construct tumour region masks... from
malignant-cell probabilities") is where spatial smoothing and thresholding
into hard regions happens, per the specification's phrasing.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.epithelial_subset import DIAGNOSTIC_RESOLUTION
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

TUMOUR_MARKER_GENES = ["EPCAM", "MET", "ERBB2", "EGFR", "KRT7"]
HPV_PROBE_PREFIX = "HPV16_"
MIN_MARKERS_PRESENT = 2


def _zscore(series: pd.Series) -> pd.Series:
    mean = series.mean(skipna=True)
    std = series.std(skipna=True)
    if series.notna().sum() == 0:
        # No data at all (e.g. every cell in this group is NaN, such as a
        # patient with no HPV panel) -- must stay NaN, not be coerced to
        # 0.0, or "no data" would silently be treated as "exactly average"
        # and wrongly pulled into the combined score's mean.
        return series
    if std == 0 or np.isnan(std):
        return pd.Series(np.where(series.notna(), 0.0, np.nan), index=series.index)
    return (series - mean) / std


def _zscore_within_patient(series: pd.Series, patient_ids: pd.Series) -> pd.Series:
    """Per-patient z-score, not a single pooled one -- required because raw
    HPV probe signal in particular varies by 2-3 orders of magnitude across
    patients for reasons unrelated to per-cell malignancy (viral copy
    number / tumour purity / section depth: median HPV16 transcript
    counts per epithelial cell range from 0 in most clinically HPV+
    patients' epithelial compartments to 38 in one, confirmed
    against raw counts). A
    single global z-score lets the one high-signal patient dominate the
    pooled mean/std, artificially pushing every other HPV+ patient's
    cells toward the negative end even though each has its own genuine
    malignant-vs-normal contrast. This is the same principle already
    applied for the same reason in `06_cell_type_annotation/03_map_external_scrna_reference.py` (per-platform standardisation
    before cross-platform comparison)
    -- the underlying question here ("which of this patient's epithelial
    cells look more malignant than that patient's own normal-epithelium
    baseline") is inherently within-patient, not a cross-patient absolute
    comparison, so standardising within patient is the correct choice for
    all four components, not only the one that first exposed the bug."""
    df = pd.DataFrame({"value": series, "patient": patient_ids})
    return df.groupby("patient", observed=True)["value"].transform(_zscore)


def compute_tumour_marker_score(
    sub: ad.AnnData, layer: str, gene_pool: list[str], rng_seed: int = get_annotation_seed()
) -> pd.Series:
    present = [g for g in TUMOUR_MARKER_GENES if g in sub.var_names]
    if len(present) < MIN_MARKERS_PRESENT:
        raise PipelineError(
            f"Only {len(present)} tumour marker gene(s) present in var_names -- below minimum {MIN_MARKERS_PRESENT}."
        )
    sc.tl.score_genes(
        sub,
        gene_list=present,
        gene_pool=gene_pool,
        layer=layer,
        score_name="_tumour_marker_score",
        random_state=rng_seed,
    )
    return sub.obs["_tumour_marker_score"]


def compute_hpv_score(
    sub: ad.AnnData,
    layer: str,
    gene_pool: list[str],
    panel_membership: pd.DataFrame,
    rng_seed: int = get_annotation_seed(),
) -> pd.Series:
    hpv_genes = [g for g in sub.var_names if g.startswith(HPV_PROBE_PREFIX)]
    result = pd.Series(np.nan, index=sub.obs_names)
    if len(hpv_genes) < MIN_MARKERS_PRESENT:
        return result

    hpv_rows = panel_membership.loc[panel_membership.index.intersection(hpv_genes)]
    sections_with_hpv = hpv_rows.columns[hpv_rows.any(axis=0)]
    hpv_cell_mask = sub.obs["section_id"].isin(sections_with_hpv).to_numpy()
    if hpv_cell_mask.sum() == 0:
        return result

    hpv_sub = sub[hpv_cell_mask].copy()
    # scanpy's score_genes bins genes by mean expression using `gene_pool`
    # alone, then looks up each `gene_list` gene's bin in that binning --
    # gene_list must therefore be a subset of gene_pool, which the
    # biological_gene-only pool (deliberately excluding HPV/CDR3 probes as
    # control genes, to avoid composition-bias distortion, `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`) is
    # not. The HPV genes being scored are added to the pool locally, for
    # this call only -- they still cannot be selected as control genes for
    # a *different* target gene elsewhere, since this extended pool is
    # never passed to any other score_genes call in this module.
    sc.tl.score_genes(
        hpv_sub,
        gene_list=hpv_genes,
        gene_pool=gene_pool + hpv_genes,
        layer=layer,
        score_name="_hpv_score",
        random_state=rng_seed,
    )
    result.loc[hpv_sub.obs_names] = hpv_sub.obs["_hpv_score"].to_numpy()
    return result


def compute_patient_clonality_score(cluster_labels: pd.Series, patient_ids: pd.Series) -> pd.Series:
    """Per-cell version of `epithelial_subset.compute_joint_cluster_patient_dominance`
    -- each cell gets its own cluster's dominant-patient fraction."""
    df = pd.DataFrame({"cluster": cluster_labels, "patient": patient_ids})
    dominance_per_cluster = df.groupby("cluster", observed=True)["patient"].apply(
        lambda s: s.value_counts(normalize=True).iloc[0]
    )
    # Leiden cluster labels are categorical dtype (scanpy's output); mapping
    # a categorical Series through a lookup can itself come back categorical
    # rather than float in pandas, which then breaks downstream numeric
    # reductions (`.mean()`) with a TypeError -- cast explicitly rather
    # than relying on `.map()`'s dtype inference.
    return df["cluster"].map(dominance_per_cluster).astype(float)


def combine_malignancy_evidence(
    tumour_marker_score: pd.Series,
    hpv_score: pd.Series,
    emt_stress_score: pd.Series,
    patient_clonality_score: pd.Series,
    patient_ids: pd.Series,
) -> pd.DataFrame:
    components = pd.DataFrame(
        {
            "tumour_marker_score": _zscore_within_patient(tumour_marker_score, patient_ids),
            "hpv_score": _zscore_within_patient(hpv_score, patient_ids),
            "emt_stress_score": _zscore_within_patient(emt_stress_score, patient_ids),
            "patient_clonality_score": _zscore_within_patient(patient_clonality_score, patient_ids),
        }
    )
    malignancy_score = components.mean(axis=1, skipna=True)
    malignancy_probability = malignancy_score.rank(pct=True)
    result = components.copy()
    result["malignancy_score"] = malignancy_score
    result["malignancy_probability"] = malignancy_probability
    return result


def build_malignancy_score_report(project_root: Path) -> dict:
    epithelial_subset_path = project_root / "data" / "objects" / "epithelial_subset.h5ad"
    panel_membership_path = (
        project_root
        / "results"
        / "tables"
        / "03_spatialdata_import"
        / "gene_panel_membership.parquet"
    )
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "data" / "derived" / "malignancy_scores.parquet"

    for p in (epithelial_subset_path, panel_membership_path, feature_annotation_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    sub = ad.read_h5ad(epithelial_subset_path)
    layer = sub.uns["primary_normalization_layer"]
    panel_membership = pd.read_parquet(panel_membership_path)
    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")
    gene_pool = feature_annotation.loc[
        feature_annotation["feature_class"] == "biological_gene", "feature_name"
    ].tolist()
    gene_pool = [g for g in gene_pool if g in sub.var_names]

    tumour_marker_score = compute_tumour_marker_score(sub, layer, gene_pool)
    hpv_score = compute_hpv_score(sub, layer, gene_pool, panel_membership)
    emt_stress_score = sub.obs[["emt_score", "stress_score"]].mean(axis=1)
    cluster_col = f"joint_leiden_res{DIAGNOSTIC_RESOLUTION}"
    patient_clonality_score = compute_patient_clonality_score(
        sub.obs[cluster_col], sub.obs["patient_id"]
    )

    result = combine_malignancy_evidence(
        tumour_marker_score,
        hpv_score,
        emt_stress_score,
        patient_clonality_score,
        sub.obs["patient_id"],
    )
    result.index = sub.obs_names
    result["patient_id"] = sub.obs["patient_id"].to_numpy()
    result["section_id"] = sub.obs["section_id"].to_numpy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    n_hpv_scored = int(hpv_score.notna().sum())
    return {
        "n_cells": len(result),
        "n_cells_hpv_scored": n_hpv_scored,
        "fraction_cells_hpv_scored": round(n_hpv_scored / len(result), 4),
        "mean_malignancy_score": round(float(result["malignancy_score"].mean()), 4),
        "malignancy_score_by_patient": {
            str(k): round(float(v), 4)
            for k, v in result.groupby("patient_id", observed=True)["malignancy_score"]
            .mean()
            .items()
        },
        "output_path": str(output_path),
    }
