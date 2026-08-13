"""Cell-level QC metric computation (`04_quality_control/00_compute_cell_level_qc_metrics.py`).

Operates on the `03_spatialdata_import/05_build_combined_analysis_object.py` combined object, which already carries
transcript_counts, control_probe_counts, control_codeword_counts,
unassigned_codeword_counts, deprecated_codeword_counts, cell_area and
nucleus_area per cell (joined from cells.parquet by the `03_spatialdata_import/01_import_each_section_to_spatialdata.py` reader)
-- this module adds detected-gene counts (from X, independently of the
pre-computed transcript_counts) and derived ratios/densities, and
cross-checks transcript_counts against the sum of X as an integrity
check rather than trusting the pre-computed field blindly.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

REQUIRED_OBS_COLUMNS = [
    "transcript_counts",
    "control_probe_counts",
    "control_codeword_counts",
    "cell_area",
    "nucleus_area",
    "section_id",
    "patient_id",
]


def compute_cell_qc_metrics(adata: ad.AnnData) -> pd.DataFrame:
    missing = [c for c in REQUIRED_OBS_COLUMNS if c not in adata.obs.columns]
    if missing:
        raise PipelineError(
            f"AnnData obs is missing required column(s) {missing} -- expected from the `03_spatialdata_import/01_import_each_section_to_spatialdata.py`, `03_spatialdata_import/05_build_combined_analysis_object.py` join."
        )

    X = adata.X
    n_genes_detected = np.asarray((X > 0).sum(axis=1)).ravel()
    counts_from_matrix = np.asarray(X.sum(axis=1)).ravel()

    metrics = adata.obs[
        [
            "section_id",
            "patient_id",
            "transcript_counts",
            "control_probe_counts",
            "control_codeword_counts",
            "cell_area",
            "nucleus_area",
        ]
    ].copy()
    metrics["n_genes_detected"] = n_genes_detected
    metrics["counts_from_expression_matrix"] = counts_from_matrix

    # 10x's cells.parquet `transcript_counts` field is gene-expression-only
    # (it does not include control_probe_counts/control_codeword_counts,
    # which are separate, sibling per-cell fields): it is exactly equal
    # to counts_from_expression_matrix for all 1,186,916 cells in this
    # project's combined object, including cells with nonzero control counts.
    # counts_discrepancy is therefore expected to be identically zero in
    # practice; it remains an integrity check (adata.X faithfully
    # reproduces cells.parquet's transcript_counts, i.e. no cells/genes lost
    # or duplicated during the `03_spatialdata_import/03_create_anndata_expression_objects.py` AnnData export), just not a check
    # on control-probe inclusion.
    metrics["counts_discrepancy"] = (
        metrics["transcript_counts"] - metrics["counts_from_expression_matrix"]
    )

    metrics["control_probe_ratio"] = metrics["control_probe_counts"] / metrics[
        "transcript_counts"
    ].replace(0, np.nan)
    metrics["control_codeword_ratio"] = metrics["control_codeword_counts"] / metrics[
        "transcript_counts"
    ].replace(0, np.nan)
    metrics["transcript_density_per_um2"] = metrics["transcript_counts"] / metrics[
        "cell_area"
    ].replace(0, np.nan)
    metrics["nucleus_to_cell_area_ratio"] = metrics["nucleus_area"] / metrics["cell_area"].replace(
        0, np.nan
    )

    return metrics


def build_cell_qc_report(combined_h5ad_path: Path, output_path: Path) -> dict:
    if not combined_h5ad_path.is_file():
        raise PipelineError(
            f"'{combined_h5ad_path}' not found. Run `03_spatialdata_import/05_build_combined_analysis_object.py` first."
        )

    adata = ad.read_h5ad(combined_h5ad_path)
    metrics = compute_cell_qc_metrics(adata)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_parquet(output_path)

    return {
        "n_cells": len(metrics),
        "n_sections": metrics["section_id"].nunique(),
        "median_transcript_counts": float(metrics["transcript_counts"].median()),
        "median_genes_detected": float(metrics["n_genes_detected"].median()),
        "median_control_probe_ratio": float(metrics["control_probe_ratio"].median()),
    }
