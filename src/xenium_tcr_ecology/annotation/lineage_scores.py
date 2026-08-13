"""Major lineage score calculation (`06_cell_type_annotation/02_score_major_lineages.py`).

Computes per-cell scores for each of `06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s `major_lineage` identities
(Epithelial_Tumour, T_cell, NK_cell, B_cell, Plasma_cell, Myeloid,
Dendritic_cell, Mast_cell, Fibroblast, Endothelial,
Perivascular_SmoothMuscle, Erythroid) using scanpy's standard `score_genes`
method, on the primary normalisation layer selected in `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`, mirroring
`05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s program-score approach exactly (same method, same exposure-gene
control pool -- the 399 `biological_gene` panel, excluding patient-specific
TCR/CDR3/HPV probes).

This step assigns SCORES only, not final cell-type labels: per the
blueprint's own caution (`06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`'s docstring) and the Cell Type Annotation completion
gate ("every cell has a label, confidence and ambiguity status"), combining
these scores with clustering, reference-mapping and spatial evidence into a
final call is explicitly `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s job, not this one.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

MIN_GENES_PER_LINEAGE = 2


def load_major_lineage_gene_sets(marker_registry: pd.DataFrame) -> dict[str, list[str]]:
    major = marker_registry[marker_registry["hierarchy_level"] == "major_lineage"]
    if len(major) == 0:
        raise PipelineError("Marker registry has no major_lineage entries.")
    gene_sets = {}
    for _, row in major.iterrows():
        genes = [g for g in str(row["markers"]).split(";") if g]
        if len(genes) < MIN_GENES_PER_LINEAGE:
            raise PipelineError(
                f"Major lineage '{row['cell_identity']}' has only {len(genes)} marker(s) -- "
                f"below the minimum of {MIN_GENES_PER_LINEAGE}."
            )
        gene_sets[row["cell_identity"]] = genes
    return gene_sets


def compute_lineage_scores(
    adata: ad.AnnData,
    layer: str,
    gene_sets: dict[str, list[str]],
    gene_pool: list[str],
    rng_seed: int = get_annotation_seed(),
) -> pd.DataFrame:
    scores = pd.DataFrame(index=adata.obs_names)
    for identity, genes in gene_sets.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < MIN_GENES_PER_LINEAGE:
            raise PipelineError(
                f"Major lineage '{identity}' has only {len(present)} marker(s) present in "
                f"adata.var_names -- below the minimum of {MIN_GENES_PER_LINEAGE}."
            )
        score_name = f"{identity}_lineage_score"
        sc.tl.score_genes(
            adata,
            gene_list=present,
            gene_pool=gene_pool,
            layer=layer,
            score_name=score_name,
            random_state=rng_seed,
        )
        scores[score_name] = adata.obs[score_name]
    return scores


def build_lineage_scores_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    marker_registry_path = project_root / "references" / "cell_type_marker_registry.tsv"
    output_path = project_root / "data" / "derived" / "lineage_scores.parquet"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )
    if not marker_registry_path.is_file():
        raise PipelineError(
            f"'{marker_registry_path}' not found. Run `06_cell_type_annotation/00_compile_marker_and_reference_registry.py` first."
        )

    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    marker_registry = pd.read_csv(marker_registry_path, sep="\t")
    gene_sets = load_major_lineage_gene_sets(marker_registry)

    gene_pool = adata.var_names[adata.var["is_exposure_gene"]].tolist()
    scores = compute_lineage_scores(adata, layer, gene_sets, gene_pool)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(output_path)

    predicted_lineage = scores.idxmax(axis=1).str.replace("_lineage_score", "", regex=False)
    lineage_counts = predicted_lineage.value_counts().to_dict()

    return {
        "n_cells": len(scores),
        "n_lineages": len(gene_sets),
        "genes_per_lineage": {k: len(v) for k, v in gene_sets.items()},
        "argmax_lineage_counts": lineage_counts,
    }
