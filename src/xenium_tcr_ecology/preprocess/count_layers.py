"""Analysis count layer construction (`05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`).

Builds `layers["counts"]` (raw, preserved explicitly), `layers["lognorm"]`
(total-count normalisation + log1p), `layers["pearson_residuals"]`
(analytic Pearson residual variance stabilisation, Lause et al. 2021), and
`layers["detected"]` (binary detection) on top of `04_quality_control/07_apply_qc_filters_with_audit_trail.py`'s QC-filtered
object, leaving `adata.X` as the untouched raw counts (Preprocessing and Normalisation's own
framing: "create analysis-ready expression representations without
removing genuine patient biology" -- multiple parallel representations, not
a single early commitment; `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` evaluates which is most appropriate
for downstream use).

Both normalisation methods compute their per-cell "exposure" (total-count
size factor for lognorm; the offset term for Pearson residuals) using ONLY
`biological_gene`-class features (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s classification), not the
full 623-feature panel. Checked against data before choosing
this over the simpler whole-panel default: the 224 patient-specific
TCR/CDR3 and HPV probes contribute a pooled mean of 6.4% of total per-cell
counts, but up to 69.6% for individual cells with strong clonal/viral
signal. Including them in the exposure computation would inflate the
apparent "size" of exactly the cells most biologically interesting to this
project (strong TCR clonotype expression), systematically suppressing their
normalised biological-gene values -- a textbook composition-bias artefact,
not a hypothetical concern. `layers["counts"]` and `layers["detected"]` are
unaffected by this choice (they are not normalised).

The Pearson residual formula here is a direct, verified generalisation of
scanpy's own `sc.experimental.pp.normalize_pearson_residuals` analytic
formula (`mu_ij = row_sum_i * col_sum_j / grand_total`), substituting a
`biological_gene`-restricted row_sum/grand_total while keeping col_sum_j
(each gene's own total, computed across all cells) unrestricted -- verified
to exactly reproduce scanpy's own output in the special case where the
"restricted" set is every gene (see
tests/unit/test_preprocess_count_layers.py).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Matches scanpy's own default (sc.experimental.pp.normalize_pearson_residuals),
# used here for direct comparability with the standard method.
PEARSON_RESIDUALS_THETA = 100.0
PEARSON_RESIDUALS_CLIP = None  # None -> scanpy's default of sqrt(n_obs)

REQUIRED_FEATURE_COLUMNS = ["feature_name", "feature_class", "in_analysis_matrix"]


def compute_size_factors(
    X, exposure_gene_mask: np.ndarray, target_sum: float | None = None
) -> tuple[np.ndarray, float]:
    """Per-cell size factors from a restricted gene subset."""
    restricted_counts = np.asarray(X[:, exposure_gene_mask].sum(axis=1)).ravel()
    if target_sum is None:
        nonzero = restricted_counts[restricted_counts > 0]
        if len(nonzero) == 0:
            raise PipelineError(
                "All cells have zero counts among the exposure gene set -- cannot compute size factors."
            )
        target_sum = float(np.median(nonzero))
    size_factors = restricted_counts / target_sum
    return size_factors, target_sum


def normalize_log1p(X, size_factors: np.ndarray):
    safe_factors = np.where(size_factors > 0, size_factors, 1.0)
    if sparse.issparse(X):
        inv = sparse.diags(1.0 / safe_factors)
        normalized = (inv @ X).tocsr()
        return normalized.log1p()
    return np.log1p(X / safe_factors[:, None])


def compute_pearson_residuals(
    X,
    exposure_gene_mask: np.ndarray,
    theta: float = PEARSON_RESIDUALS_THETA,
    clip: float | None = PEARSON_RESIDUALS_CLIP,
) -> np.ndarray:
    if theta <= 0:
        raise PipelineError("Pearson residuals require theta > 0.")

    X_dense = X.toarray() if sparse.issparse(X) else np.asarray(X)
    n_obs = X_dense.shape[0]
    if clip is None:
        clip = np.sqrt(n_obs)

    restricted_row_sums = X_dense[:, exposure_gene_mask].sum(axis=1, keepdims=True)
    restricted_grand_total = restricted_row_sums.sum()
    col_sums = X_dense.sum(axis=0, keepdims=True)

    mu = restricted_row_sums @ col_sums / restricted_grand_total
    residuals = (X_dense - mu) / np.sqrt(mu + mu**2 / theta)
    return np.clip(residuals, -clip, clip)


def build_analysis_count_layers(adata: ad.AnnData, feature_annotation: pd.DataFrame) -> ad.AnnData:
    missing = [c for c in REQUIRED_FEATURE_COLUMNS if c not in feature_annotation.columns]
    if missing:
        raise PipelineError(f"feature_annotation is missing required column(s) {missing}.")

    bio_gene_names = set(
        feature_annotation.loc[
            feature_annotation["feature_class"] == "biological_gene", "feature_name"
        ]
    )
    exposure_mask = np.array([g in bio_gene_names for g in adata.var_names])
    if exposure_mask.sum() == 0:
        raise PipelineError(
            "No biological_gene features found in adata.var_names -- cannot compute an exposure basis."
        )

    adata.layers["counts"] = adata.X.copy()

    size_factors, target_sum = compute_size_factors(adata.X, exposure_mask)
    adata.obs["size_factor"] = size_factors
    adata.uns["lognorm_target_sum"] = target_sum
    adata.layers["lognorm"] = normalize_log1p(adata.X, size_factors)

    adata.layers["pearson_residuals"] = compute_pearson_residuals(adata.X, exposure_mask)
    adata.uns["pearson_residuals_theta"] = PEARSON_RESIDUALS_THETA

    detected = adata.X.copy()
    if sparse.issparse(detected):
        detected = (detected > 0).astype(np.float32)
    else:
        detected = (np.asarray(detected) > 0).astype(np.float32)
    adata.layers["detected"] = detected

    adata.var["is_exposure_gene"] = exposure_mask
    return adata


def build_count_layers_report(project_root: Path) -> dict:
    qc_filtered_path = project_root / "data" / "objects" / "qc_filtered.h5ad"
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "data" / "objects" / "analysis_ready.h5ad"

    if not qc_filtered_path.is_file():
        raise PipelineError(
            f"'{qc_filtered_path}' not found. Run `04_quality_control/07_apply_qc_filters_with_audit_trail.py` first."
        )
    if not feature_annotation_path.is_file():
        raise PipelineError(
            f"'{feature_annotation_path}' not found. Run `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py` first."
        )

    adata = ad.read_h5ad(qc_filtered_path)
    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")

    adata = build_analysis_count_layers(adata, feature_annotation)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_exposure_genes": int(adata.var["is_exposure_gene"].sum()),
        "layers": sorted(adata.layers.keys()),
        "lognorm_target_sum": float(adata.uns["lognorm_target_sum"]),
        "median_size_factor": float(np.median(adata.obs["size_factor"])),
    }
