"""Normalisation strategy evaluation (`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`).

Compares `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`'s four layers (counts, lognorm, pearson_residuals,
detected) on two axes named explicitly in the blueprint: replicate
stability and negative-control-probe behaviour.

Runs in Python, not R, despite the blueprint assigning this phase step an
`.R` script: the comparison needs direct access to `analysis_ready.h5ad`
(1.12M cells x 623 genes x 4 layers), and no R HDF5/AnnData reader is
available in this project's environment, the same R-side interop
constraint already hit in `04_quality_control/08_assess_replicate_concordance.R`. This module is invoked by
`scripts/05_preprocessing_and_normalisation/_02_compute_normalization_benchmark_metrics.py`
(a helper, not its own numbered blueprint phase step) via a subprocess call
from `02_evaluate_normalisation_strategies.R`, which owns the actual
blueprint-mandated PDF report -- the same "Python computes, R reports"
split already established for the `03_spatialdata_import/06_export_r_interoperability_objects.py` R-export bridge.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

from xenium_tcr_ecology.infra.exceptions import PipelineError

METHODS = ["counts", "lognorm", "pearson_residuals", "detected"]

# Human decision, already made and signed off
# (log-normalisation selected as the primary layer, based on the
# documented evidence comparison there). This module reports the
# comparison; it does not re-decide
# the winner each run. `apply_primary_normalization_layer_decision` below
# is what actually persists that decision onto analysis_ready.h5ad, for
# `05_preprocessing_and_normalisation/03_calculate_program_scores.py`
# onward to read.
PRIMARY_NORMALIZATION_LAYER = "lognorm"


def _layer_as_dense(adata: ad.AnnData, layer: str, gene_mask: np.ndarray) -> np.ndarray:
    values = adata.layers[layer][:, gene_mask]
    return values.toarray() if sparse.issparse(values) else np.asarray(values)


def _pseudobulk_profile(adata: ad.AnnData, layer: str, gene_mask: np.ndarray) -> np.ndarray:
    values = _layer_as_dense(adata, layer, gene_mask)
    if layer == "counts":
        return np.log1p(values.sum(axis=0))
    return values.mean(axis=0)


def compute_replicate_stability(
    analysis_ready_path: Path, exclusion_log_path: Path
) -> pd.DataFrame:
    if not analysis_ready_path.is_file():
        raise PipelineError(
            f"'{analysis_ready_path}' not found. Run `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py` first."
        )
    if not exclusion_log_path.is_file():
        raise PipelineError(
            f"'{exclusion_log_path}' not found. Run `04_quality_control/07_apply_qc_filters_with_audit_trail.py` first."
        )

    adata = ad.read_h5ad(analysis_ready_path)
    exclusion_log = pd.read_csv(exclusion_log_path, sep="\t")
    exclusion_log = exclusion_log.set_index("cell_id")

    section_to_patient = adata.obs.groupby("section_id", observed=True)["patient_id"].first()
    replicate_sections = adata.obs.loc[adata.obs["is_technical_replicate"], "section_id"].unique()
    patients = adata.obs.loc[
        adata.obs["section_id"].isin(replicate_sections), "patient_id"
    ].unique()

    gene_mask = adata.var["is_exposure_gene"].to_numpy()

    rows = []
    for patient_id in sorted(patients):
        sections = sorted(section_to_patient.index[section_to_patient == patient_id])
        if len(sections) != 2:
            raise PipelineError(
                f"Patient '{patient_id}' has {len(sections)} replicate section(s), not 2."
            )
        section1, section2 = sections

        mask1 = (adata.obs["section_id"] == section1) & (
            exclusion_log.reindex(adata.obs_names)["qc_pass"].to_numpy()
        )
        mask2 = (adata.obs["section_id"] == section2) & (
            exclusion_log.reindex(adata.obs_names)["qc_pass"].to_numpy()
        )

        adata1 = adata[mask1.to_numpy()]
        adata2 = adata[mask2.to_numpy()]

        row = {"patient_id": patient_id, "section1": section1, "section2": section2}
        for method in METHODS:
            profile1 = _pseudobulk_profile(adata1, method, gene_mask)
            profile2 = _pseudobulk_profile(adata2, method, gene_mask)
            r, _ = stats.pearsonr(profile1, profile2)
            row[f"{method}_replicate_r"] = r
        rows.append(row)

    return pd.DataFrame(rows)


def compute_technical_noise_correlation(
    analysis_ready_path: Path, cell_qc_metrics_path: Path
) -> pd.DataFrame:
    if not analysis_ready_path.is_file():
        raise PipelineError(
            f"'{analysis_ready_path}' not found. Run `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py` first."
        )
    if not cell_qc_metrics_path.is_file():
        raise PipelineError(
            f"'{cell_qc_metrics_path}' not found. Run `04_quality_control/00_compute_cell_level_qc_metrics.py` first."
        )

    adata = ad.read_h5ad(analysis_ready_path)
    cell_qc = pd.read_parquet(cell_qc_metrics_path)
    control_probe_ratio = (
        cell_qc.reindex(adata.obs_names)["control_probe_ratio"].fillna(0.0).to_numpy()
    )

    gene_mask = adata.var["is_exposure_gene"].to_numpy()

    rows = []
    for method in METHODS:
        values = _layer_as_dense(adata, method, gene_mask)
        per_cell_mean = values.mean(axis=1)
        rho, p_value = stats.spearmanr(per_cell_mean, control_probe_ratio)
        rows.append(
            {"method": method, "spearman_rho_vs_control_probe_ratio": rho, "p_value": p_value}
        )

    return pd.DataFrame(rows)


def apply_primary_normalization_layer_decision(analysis_ready_path: Path) -> bool:
    """Persists the already-governance-decided primary normalisation layer
    (see `PRIMARY_NORMALIZATION_LAYER` above) onto analysis_ready.h5ad's
    `uns`, so `03_calculate_program_scores.py` and
    `preprocess/release_freeze.py` can read it. Idempotent: does nothing
    (returns False) if already set, so repeat runs don't rewrite the file
    unnecessarily. Returns True if it wrote the file.

    Found missing during the second Vast.ai clean-room run: nothing in this
    project previously set this attribute -- `02_evaluate_normalisation_strategies.R`
    only reports the comparison (see its own docstring and this module's),
    and no R HDF5/AnnData writer is available to it either.
    """
    adata = ad.read_h5ad(analysis_ready_path)
    if adata.uns.get("primary_normalization_layer") == PRIMARY_NORMALIZATION_LAYER:
        return False
    adata.uns["primary_normalization_layer"] = PRIMARY_NORMALIZATION_LAYER
    adata.write_h5ad(analysis_ready_path)
    return True


def build_normalization_benchmark_summary(project_root: Path, output_path: Path) -> dict:
    analysis_ready_path = project_root / "data" / "objects" / "analysis_ready.h5ad"
    exclusion_log_path = project_root / "data" / "derived" / "exclusion_log.tsv"
    cell_qc_metrics_path = project_root / "data" / "derived" / "cell_qc_metrics.parquet"

    replicate_stability = compute_replicate_stability(analysis_ready_path, exclusion_log_path)
    technical_noise = compute_technical_noise_correlation(analysis_ready_path, cell_qc_metrics_path)
    apply_primary_normalization_layer_decision(analysis_ready_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    replicate_stability.to_parquet(
        output_path.with_name(output_path.stem + "_replicate_stability.parquet")
    )
    technical_noise.to_parquet(output_path.with_name(output_path.stem + "_technical_noise.parquet"))

    median_replicate_r = {
        method: float(replicate_stability[f"{method}_replicate_r"].median()) for method in METHODS
    }
    abs_technical_noise_rho = dict(
        zip(technical_noise["method"], technical_noise["spearman_rho_vs_control_probe_ratio"].abs())
    )

    return {
        "n_replicate_pairs": len(replicate_stability),
        "median_replicate_r_by_method": median_replicate_r,
        "abs_technical_noise_rho_by_method": {
            k: float(v) for k, v in abs_technical_noise_rho.items()
        },
    }
