"""Epithelial cell subsetting and reclustering (`07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py`).

Characterises epithelial heterogeneity separately from immune/stromal
structure: the 397,247 cells `06_cell_type_annotation/06_integrate_annotation_evidence.py` called `final_lineage ==
"Epithelial_Tumour"` are subset out of the primary analysis matrix and
reclustered on their own, at both patient (`within_patient_leiden`) and
joint (`joint_leiden_res*`) scope, reusing `06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`'s own generic
clustering functions (`compute_joint_clustering`,
`compute_within_patient_clustering`) rather than duplicating them --
nothing about those functions is specific to the full dataset, only to
"some AnnData with a `patient_id` obs column and a gene-expression layer."

Clustering genes are restricted to the same `biological_gene` set Phase
6.01 already used (excludes the 216 CDR3 and 8 HPV probe features), for
the same documented reason:
these are extremely sparse and panel-variable by patient, and would let
unsupervised clustering trivially recover patient/HPV-status identity
rather than discover genuine transcriptional structure -- exactly the
confound this subsetting step needs to avoid, doubly so here since
malignant-cell clustering is already expected to be patient-dominated for
biological reasons (see below), and an artificial, panel-driven
patient signal would be indistinguishable from the genuine one.

**Joint vs within-patient clustering is expected to diverge sharply here,
and this is itself the correct, literature-consistent finding, not an
artefact to fix:** malignant/tumour epithelial cells are well established
in the HNSCC single-cell literature (Puram et al. 2017, Cell -- the
foundational HNSCC single-cell malignant-heterogeneity study) to cluster
predominantly by patient of origin (driven by patient-specific
copy-number/genotype), not by a shared cross-patient malignant transcriptional
state, in sharp contrast to immune/stromal populations which do cluster by
cell type across patients (already observed for the *whole* dataset in
`06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`, `06_cell_type_annotation/06_integrate_annotation_evidence.py`). `build_epithelial_subset_report` therefore reports an
explicit joint-cluster/patient-dominance diagnostic (the mean fraction of
each joint cluster's cells belonging to that cluster's single most common
patient) precisely so this expected pattern is measured and recorded here,
not silently assumed.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from xenium_tcr_ecology.annotation.clustering import (
    compute_joint_clustering,
    compute_within_patient_clustering,
    select_clustering_gene_mask,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError

EPITHELIAL_LINEAGE_LABEL = "Epithelial_Tumour"
# Matches `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s own choice of a moderate, non-extreme resolution
# (integrate_evidence.py's CLUSTER_RESOLUTION_COL) -- reused here as the
# single resolution reported for the patient-dominance diagnostic below.
DIAGNOSTIC_RESOLUTION = 0.5


def subset_epithelial_cells(adata: ad.AnnData, final_annotations: pd.DataFrame) -> ad.AnnData:
    if EPITHELIAL_LINEAGE_LABEL not in final_annotations["final_lineage"].unique():
        raise PipelineError(f"No cells with final_lineage == '{EPITHELIAL_LINEAGE_LABEL}' found.")

    epithelial_ids = final_annotations.index[
        final_annotations["final_lineage"] == EPITHELIAL_LINEAGE_LABEL
    ]
    common_ids = adata.obs_names.intersection(epithelial_ids)
    if len(common_ids) == 0:
        raise PipelineError("No overlap between adata.obs_names and epithelial cell IDs.")

    sub = adata[common_ids].copy()
    sub.obs["confidence"] = final_annotations.loc[common_ids, "confidence"]
    sub.obs["is_ambiguous"] = final_annotations.loc[common_ids, "is_ambiguous"]
    return sub


def compute_joint_cluster_patient_dominance(
    cluster_labels: pd.Series, patient_ids: pd.Series
) -> float:
    """Mean, over joint clusters, of the fraction of a cluster's cells
    belonging to that cluster's single most common patient -- a direct,
    interpretable measure of how patient-dominated the joint clustering
    structure is (1.0 = every cluster is a single patient; low values =
    clusters mix patients freely). This cluster-unweighted mean treats a
    37-cell cluster and a 51,000-cell cluster identically, so pair it with
    `compute_cell_weighted_patient_dominance` below when cluster sizes vary
    a lot (checked on the data: they do here, 37-51,223 cells,
    and the two metrics diverge meaningfully as a resultTumour Epithelium Characterisation companion entry).
    """
    df = pd.DataFrame({"cluster": cluster_labels, "patient": patient_ids})
    dominance_per_cluster = df.groupby("cluster", observed=True)["patient"].apply(
        lambda s: s.value_counts(normalize=True).iloc[0]
    )
    return float(dominance_per_cluster.mean())


def compute_cell_weighted_patient_dominance(
    cluster_labels: pd.Series, patient_ids: pd.Series
) -> float:
    """The same per-cluster dominance fraction as
    `compute_joint_cluster_patient_dominance`, but averaged weighted by
    cluster size (i.e. equivalently: the fraction of ALL cells that sit in
    their own cluster's single most common patient) -- answers "how
    patient-dominated is the clustering structure that most cells actually
    experience," which the unweighted per-cluster mean does not when
    cluster sizes are highly unequal."""
    df = pd.DataFrame({"cluster": cluster_labels, "patient": patient_ids})
    dominance_per_cluster = df.groupby("cluster", observed=True)["patient"].apply(
        lambda s: s.value_counts(normalize=True).iloc[0]
    )
    cluster_sizes = df.groupby("cluster", observed=True).size()
    return float((dominance_per_cluster * cluster_sizes).sum() / cluster_sizes.sum())


def build_epithelial_subset_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "data" / "objects" / "epithelial_subset.h5ad"

    for p in (matrix_path, final_annotations_path, feature_annotation_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    adata = ad.read_h5ad(matrix_path)
    final_annotations = pd.read_parquet(final_annotations_path)
    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")

    sub = subset_epithelial_cells(adata, final_annotations)
    layer = sub.uns["primary_normalization_layer"]
    gene_mask = select_clustering_gene_mask(sub.var_names, feature_annotation)

    joint = compute_joint_clustering(sub, layer, gene_mask)
    within_patient = compute_within_patient_clustering(sub, layer, gene_mask)

    for col in joint.columns:
        sub.obs[col] = joint[col].to_numpy()
    sub.obs["within_patient_leiden"] = within_patient.to_numpy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sub.write_h5ad(output_path)

    diagnostic_col = f"joint_leiden_res{DIAGNOSTIC_RESOLUTION}"
    patient_dominance = compute_joint_cluster_patient_dominance(
        sub.obs[diagnostic_col], sub.obs["patient_id"]
    )
    cell_weighted_dominance = compute_cell_weighted_patient_dominance(
        sub.obs[diagnostic_col], sub.obs["patient_id"]
    )

    n_clusters_by_resolution = {col: int(sub.obs[col].nunique()) for col in joint.columns}
    return {
        "n_epithelial_cells": sub.n_obs,
        "n_patients_represented": int(sub.obs["patient_id"].nunique()),
        "fraction_ambiguous": round(float(sub.obs["is_ambiguous"].mean()), 4),
        "n_clustering_genes": int(gene_mask.sum()),
        "n_clusters_by_resolution": n_clusters_by_resolution,
        "n_within_patient_clusters": int(sub.obs["within_patient_leiden"].nunique()),
        "n_cells_with_within_patient_label": int(sub.obs["within_patient_leiden"].notna().sum()),
        "joint_cluster_patient_dominance": round(patient_dominance, 4),
        "joint_cluster_patient_dominance_cell_weighted": round(cell_weighted_dominance, 4),
        "diagnostic_resolution": DIAGNOSTIC_RESOLUTION,
        "output_path": str(output_path),
    }
