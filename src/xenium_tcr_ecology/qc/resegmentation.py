"""Independent alternative resegmentation of a representative section subset
(`04_quality_control/05_resegment_reference_subset.py`).

Purpose, per the specification: test whether headline spatial results are
robust to segmentation choice, not produce a second "better" segmentation.
This module therefore deliberately implements a different
algorithm class from the primary segmentation, not a variant of it.

**Primary segmentation** (already in use throughout this project, Phase
3.01 onward): 10x's own multimodal cell-boundary algorithm, which uses
actual membrane/interior protein and RNA staining to draw an irregular
polygon per cell, expanded outward from the nucleus boundary (see
`segmentation_quality.py`'s module docstring). This is morphology-aware.

**Alternative method (this module):** morphology-free nearest-nucleus-
centroid transcript reassignment -- each transcript is assigned to
whichever nucleus centroid is closest, capped at `MAX_EXPANSION_RADIUS_UM`.
This is a standard, field-precedented baseline (it is 10x's own
documented alternative "nucleus expansion" Xenium onboard preset, used
when boundary/interior staining is unavailable or as a robustness
comparator), implementable entirely with this project's existing
dependencies (`scipy.spatial.cKDTree`, no new heavy tool such as Baysor),
and its radius is a vendor-documented default (15um) rather than an
arbitrary choice.

**Representative section subset:** a genuine judgment call, made and
documented here rather
than deferred: the sparsest, median, and densest of the 18 sections by
`04_quality_control/00_compute_cell_level_qc_metrics.py` cell count, each from a different patient
(`REPRESENTATIVE_SECTIONS`) -- chosen because tissue density is the axis
most relevant to segmentation robustness (denser tissue means more
segmentation ambiguity/crowding, the harder case for this exact kind of
comparison), and picking sections from three different patients avoids the
subset being dominated by one patient's tissue characteristics.

**Feature scope:** restricted to the 623-feature "combined panel" already
used for the primary analysis matrix (`biological_gene` + `tcr_cdr3_probe`
+ `hpv_probe`, `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s `FEATURE_CLASSES_IN_ANALYSIS_MATRIX`) --
negative-control/unassigned-codeword transcripts are excluded from both
sides of the comparison so the concordance metrics below are not diluted
by features the primary matrix itself never counted.

**QV filtering:** deliberately not applied here either, for direct
consistency with the `04_quality_control/07_apply_qc_filters_with_audit_trail.py` decision not to QV-filter the primary matrix -- this keeps the
comparison isolated to the one variable actually under test (the
segmentation/reassignment rule), not compounded with an unrelated
filtering-policy difference.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as sd
from scipy.spatial import cKDTree

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.feature_classification import (
    FEATURE_CLASSES_IN_ANALYSIS_MATRIX,
    classify_feature,
)

# Data-selected subset spanning the observed cell-density range (see
# module docstring):
# sparsest (P19_run1, 23,496 cells), median (P12_run2, 43,759 cells),
# densest (P09_run2, 268,328 cells) of the 18 sections by `04_quality_control/00_compute_cell_level_qc_metrics.py`
# cell count, each from a different patient.
REPRESENTATIVE_SECTIONS = ["P19_run1", "P12_run2", "P09_run2"]

# 10x Genomics' own documented default expansion distance for Xenium's
# alternative "nucleus expansion" onboard segmentation preset (see module
# docstring).
MAX_EXPANSION_RADIUS_UM = 15.0

UNASSIGNED_CELL_ID_SENTINEL = "UNASSIGNED"


def _is_biological_feature(name: str) -> bool:
    return classify_feature(name) in FEATURE_CLASSES_IN_ANALYSIS_MATRIX


def reassign_transcripts_to_nearest_nucleus(
    transcript_xy: np.ndarray,
    nucleus_ids: np.ndarray,
    nucleus_centroids_xy: np.ndarray,
    max_radius_um: float = MAX_EXPANSION_RADIUS_UM,
) -> np.ndarray:
    """Pure, testable nearest-neighbour reassignment -- factored out so the
    core logic is testable with plain coordinate arrays, not only via a
    full SpatialData store. Returns an object array the same length as
    `transcript_xy`: the assigned nucleus_id, or None for transcripts
    farther than `max_radius_um` from every nucleus (background)."""
    if len(nucleus_centroids_xy) == 0:
        return np.full(len(transcript_xy), None, dtype=object)
    if len(transcript_xy) == 0:
        return np.array([], dtype=object)
    tree = cKDTree(nucleus_centroids_xy)
    dist, idx = tree.query(transcript_xy, k=1)
    return np.where(dist <= max_radius_um, nucleus_ids[idx], None)


def resegment_section(zarr_path: Path) -> tuple[ad.AnnData, pd.DataFrame]:
    """Returns (resegmented nucleus x gene AnnData, per-transcript
    primary-vs-alternative assignment table for concordance analysis)."""
    if not zarr_path.exists():
        raise PipelineError(
            f"'{zarr_path}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    sdata = sd.read_zarr(zarr_path)
    nucleus_shapes = sdata["nucleus_boundaries"]
    if len(nucleus_shapes) == 0:
        raise PipelineError(f"'{zarr_path}': no nucleus boundaries found.")

    nucleus_ids = nucleus_shapes.index.to_numpy()
    centroids = nucleus_shapes.geometry.centroid
    nucleus_xy = np.column_stack([centroids.x.to_numpy(), centroids.y.to_numpy()])

    pts = sdata["transcripts"][["x", "y", "feature_name", "cell_id"]].compute()
    pts = pts[pts["feature_name"].map(_is_biological_feature)].copy()

    transcript_xy = pts[["x", "y"]].to_numpy(dtype=np.float64)
    pts["reassigned_nucleus_id"] = reassign_transcripts_to_nearest_nucleus(
        transcript_xy, nucleus_ids, nucleus_xy
    )

    assigned = pts[pts["reassigned_nucleus_id"].notna()]
    counts = (
        assigned.groupby(["reassigned_nucleus_id", "feature_name"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    matrix = counts.pivot(index="reassigned_nucleus_id", columns="feature_name", values="n")
    # Every nucleus is represented, including those that received zero
    # reassigned transcripts under the alternative method.
    matrix = matrix.reindex(index=nucleus_ids, columns=matrix.columns).fillna(0.0)

    adata = ad.AnnData(
        X=matrix.to_numpy(dtype=np.float32),
        obs=pd.DataFrame(index=matrix.index.astype(str)),
        var=pd.DataFrame(index=matrix.columns.astype(str)),
    )
    adata.obs["resegmented_total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()

    return adata, pts[["cell_id", "reassigned_nucleus_id"]]


def summarize_reassignment_concordance(assignment_table: pd.DataFrame) -> dict:
    """Direct per-transcript comparison of the primary segmentation's own
    `cell_id` assignment against this module's `reassigned_nucleus_id` --
    both are in the same section-local ID space (confirmed: cell_boundaries
    and nucleus_boundaries share an index convention, `04_quality_control/03_assess_segmentation_quality.py`), so this
    is a literal, not approximate, per-transcript agreement check."""
    primary = assignment_table["cell_id"]
    alt = assignment_table["reassigned_nucleus_id"]
    primary_assigned = primary != UNASSIGNED_CELL_ID_SENTINEL
    alt_assigned = alt.notna()

    both_assigned = primary_assigned & alt_assigned
    n_both_assigned = int(both_assigned.sum())
    n_same_cell = int((both_assigned & (primary == alt)).sum())

    n = len(assignment_table)
    return {
        "n_transcripts": n,
        "fraction_concordant_same_cell": (
            round(n_same_cell / n_both_assigned, 4) if n_both_assigned else None
        ),
        "fraction_primary_assigned_alt_background": round(
            float((primary_assigned & ~alt_assigned).mean()), 4
        ),
        "fraction_primary_unassigned_alt_recovered": round(
            float((~primary_assigned & alt_assigned).mean()), 4
        ),
        "fraction_both_unassigned": round(float((~primary_assigned & ~alt_assigned).mean()), 4),
    }


def build_resegmentation_report(project_root: Path) -> dict:
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    output_dir = project_root / "data" / "objects" / "resegmented_subset"
    summary_path = project_root / "reports" / "qc" / "resegmentation_concordance.tsv"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )

    import anndata as ad_io

    primary_adata = ad_io.read_h5ad(matrix_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    section_rows = []
    for section_id in REPRESENTATIVE_SECTIONS:
        zarr_path = spatialdata_root / f"{section_id}.zarr"
        reseg_adata, assignment_table = resegment_section(zarr_path)
        reseg_adata.write_h5ad(output_dir / f"{section_id}.h5ad")

        concordance = summarize_reassignment_concordance(assignment_table)

        global_ids = [f"{section_id}_{local_id}" for local_id in reseg_adata.obs_names]
        reseg_totals = pd.Series(
            reseg_adata.obs["resegmented_total_counts"].to_numpy(), index=global_ids
        )

        # AnnData subsetting (not `.obs[mask]`, which returns a bare
        # DataFrame with no `.X`) so raw counts stay reachable below.
        section_primary = primary_adata[primary_adata.obs["section_id"] == section_id]
        common_ids = reseg_totals.index.intersection(section_primary.obs_names)
        primary_totals = np.asarray(primary_adata[common_ids].layers["counts"].sum(axis=1)).ravel()
        total_count_correlation = (
            float(np.corrcoef(primary_totals, reseg_totals.loc[common_ids].to_numpy())[0, 1])
            if len(common_ids) > 1
            else None
        )

        shared_genes = reseg_adata.var_names.intersection(primary_adata.var_names)
        primary_pseudobulk = np.asarray(
            section_primary[:, shared_genes].layers["counts"].sum(axis=0)
        ).ravel()
        reseg_pseudobulk = np.asarray(reseg_adata[:, shared_genes].X.sum(axis=0)).ravel()
        pseudobulk_correlation = (
            float(np.corrcoef(primary_pseudobulk, reseg_pseudobulk)[0, 1])
            if len(shared_genes) > 1
            else None
        )

        section_rows.append(
            {
                "section_id": section_id,
                "n_nuclei": reseg_adata.n_obs,
                "n_common_cells_for_correlation": len(common_ids),
                "n_shared_genes": len(shared_genes),
                "total_count_correlation": (
                    round(total_count_correlation, 4)
                    if total_count_correlation is not None
                    else None
                ),
                "pseudobulk_correlation": (
                    round(pseudobulk_correlation, 4) if pseudobulk_correlation is not None else None
                ),
                **concordance,
            }
        )

    summary_df = pd.DataFrame(section_rows)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(summary_path, sep="\t", index=False)

    return {
        "sections_processed": len(section_rows),
        "output_dir": str(output_dir),
        "summary_path": str(summary_path),
        "median_total_count_correlation": round(
            float(summary_df["total_count_correlation"].median()), 4
        ),
        "median_pseudobulk_correlation": round(
            float(summary_df["pseudobulk_correlation"].median()), 4
        ),
        "median_fraction_concordant_same_cell": round(
            float(summary_df["fraction_concordant_same_cell"].median()), 4
        ),
    }
