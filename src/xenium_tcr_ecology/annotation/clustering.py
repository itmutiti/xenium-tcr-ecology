"""Multi-resolution clustering, within-patient and jointly (`06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`).

Generates exploratory clustering structure at several Leiden resolutions
on the pooled ("jointly") dataset, plus a single-resolution clustering
computed independently within each patient -- a diagnostic for whether the
joint clustering's structure is consistent with each patient's own local
structure, or is instead dominated by patient-specific batch effects.

Clustering is restricted to `biological_gene` features (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s
classification), explicitly excluding the 224 patient-specific TCR/CDR3
and HPV probes: these are extremely sparse, present in only a subset of
patients by panel design, and would tend to dominate clustering by
clonotype/HPV-status rather than by cell type if included -- the same
exposure-gene-set reasoning already applied in `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`'s normalisation
and `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s program-score control-gene pool.

Cluster labels are exploratory structure, not cell-type calls: nothing
here is treated as a final identity (that is `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s job, integrating
this with marker/reference/spatial evidence) -- matching the blueprint's
own explicit caution against conflating clusters with cell types.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

JOINT_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.2]
WITHIN_PATIENT_RESOLUTION = 0.8
N_PCS = 30
N_NEIGHBORS = 15
RNG_SEED = get_annotation_seed()
HVG_N_TOP_GENES = 75
MIN_CELLS_FOR_WITHIN_PATIENT_CLUSTERING = 500


def select_clustering_gene_mask(
    var_names: pd.Index, feature_annotation: pd.DataFrame
) -> np.ndarray:
    bio_genes = set(
        feature_annotation.loc[
            feature_annotation["feature_class"] == "biological_gene", "feature_name"
        ]
    )
    mask = np.array([g in bio_genes for g in var_names])
    if mask.sum() == 0:
        raise PipelineError(
            "No biological_gene features found in var_names -- cannot select clustering genes."
        )
    return mask


def _prepare_for_clustering(
    adata: ad.AnnData, layer: str, gene_mask: np.ndarray, n_top_genes: int
) -> ad.AnnData:
    sub = adata[:, gene_mask].copy()
    sub.X = sub.layers[layer]
    n_top = min(n_top_genes, sub.n_vars)
    sc.pp.highly_variable_genes(sub, n_top_genes=n_top)
    sub = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=10)
    return sub


def compute_joint_clustering(
    adata: ad.AnnData,
    layer: str,
    gene_mask: np.ndarray,
    resolutions: list[float] = JOINT_RESOLUTIONS,
    n_pcs: int = N_PCS,
    n_neighbors: int = N_NEIGHBORS,
    n_top_genes: int = HVG_N_TOP_GENES,
    rng_seed: int = RNG_SEED,
) -> pd.DataFrame:
    sub = _prepare_for_clustering(adata, layer, gene_mask, n_top_genes)
    sc.tl.pca(sub, n_comps=min(n_pcs, sub.n_vars - 1), svd_solver="arpack", random_state=rng_seed)
    sc.pp.neighbors(sub, n_neighbors=n_neighbors, random_state=rng_seed)

    results = pd.DataFrame(index=adata.obs_names)
    for res in resolutions:
        key = f"joint_leiden_res{res}"
        sc.tl.leiden(
            sub,
            resolution=res,
            key_added=key,
            flavor="igraph",
            n_iterations=2,
            random_state=rng_seed,
        )
        results[key] = sub.obs[key].reindex(adata.obs_names).to_numpy()
    return results


def compute_within_patient_clustering(
    adata: ad.AnnData,
    layer: str,
    gene_mask: np.ndarray,
    resolution: float = WITHIN_PATIENT_RESOLUTION,
    n_pcs: int = N_PCS,
    n_neighbors: int = N_NEIGHBORS,
    n_top_genes: int = HVG_N_TOP_GENES,
    rng_seed: int = RNG_SEED,
    min_cells: int = MIN_CELLS_FOR_WITHIN_PATIENT_CLUSTERING,
) -> pd.Series:
    result = pd.Series(index=adata.obs_names, dtype=object)
    for patient_id in sorted(adata.obs["patient_id"].unique()):
        patient_mask = (adata.obs["patient_id"] == patient_id).to_numpy()
        n_cells = int(patient_mask.sum())
        if n_cells < min_cells:
            continue
        sub_adata = adata[patient_mask]
        sub = _prepare_for_clustering(sub_adata, layer, gene_mask, n_top_genes)
        sc.tl.pca(
            sub,
            n_comps=min(n_pcs, sub.n_vars - 1, sub.n_obs - 1),
            svd_solver="arpack",
            random_state=rng_seed,
        )
        sc.pp.neighbors(sub, n_neighbors=min(n_neighbors, sub.n_obs - 1), random_state=rng_seed)
        sc.tl.leiden(
            sub,
            resolution=resolution,
            key_added="wp_leiden",
            flavor="igraph",
            n_iterations=2,
            random_state=rng_seed,
        )
        result.loc[sub_adata.obs_names] = [f"{patient_id}_{c}" for c in sub.obs["wp_leiden"]]
    return result


def build_clustering_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "data" / "derived" / "clustering_assignments.parquet"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )
    if not feature_annotation_path.is_file():
        raise PipelineError(
            f"'{feature_annotation_path}' not found. Run `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py` first."
        )

    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")
    gene_mask = select_clustering_gene_mask(adata.var_names, feature_annotation)

    joint = compute_joint_clustering(adata, layer, gene_mask)
    within_patient = compute_within_patient_clustering(adata, layer, gene_mask)

    result = joint.copy()
    result["within_patient_leiden"] = within_patient

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    n_clusters_by_resolution = {
        col: int(result[col].nunique())
        for col in result.columns
        if col.startswith("joint_leiden_res")
    }
    return {
        "n_cells": len(result),
        "n_clustering_genes": int(gene_mask.sum()),
        "joint_resolutions": JOINT_RESOLUTIONS,
        "n_clusters_by_resolution": n_clusters_by_resolution,
        "n_within_patient_clusters": int(result["within_patient_leiden"].nunique()),
        "n_cells_with_within_patient_label": int(result["within_patient_leiden"].notna().sum()),
    }
