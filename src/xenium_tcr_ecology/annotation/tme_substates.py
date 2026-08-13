"""Myeloid/stromal substate input preparation (`06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R` helper).

Computes marker-score evidence for 5 compartments' substates -- Myeloid
(Macrophage vs Monocyte), Dendritic_cell (cDC vs pDC vs Mature_DC),
Fibroblast (Activated_CAF vs Resting_Fibroblast), Endothelial (Blood vs
Lymphatic), Perivascular_SmoothMuscle (Pericyte vs Smooth_muscle) -- all
drawn from genes already validated as present and lineage-appropriate in
`06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s marker registry, not a fresh, unverified curation.

Unlike `06_cell_type_annotation/04_resolve_t_cell_substates.R` (T-cell substates), there is no external reference
dataset for this compartment: GSE287301 (`06_cell_type_annotation/03_map_external_scrna_reference.py`) is T cells only. This
module is therefore panel-marker-only, with no reference-transfer
cross-check available -- a documented scope difference from 6.04, not
an oversight.

Writes a plain parquet for 05_resolve_myeloid_and_stromal_substates.R to
consume directly (R cannot read analysis_ready.h5ad -- the same "Python
computes, R decides/reports" split as `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`, `06_cell_type_annotation/04_resolve_t_cell_substates.R`).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

# Every marker below is a subset of `06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s already-validated
# major_lineage marker sets (references/cell_type_marker_registry.tsv),
# split into finer, panel-supported subtypes -- not independently curated
# from scratch. ACTA2 is deliberately reused from Perivascular_SmoothMuscle
# for Activated_CAF: `06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s registry already documents this overlap
# explicitly ("ACTA2 and PDGFRB are also myofibroblast/activated-stroma
# markers").
SUBSTATE_MARKERS: dict[str, dict[str, list[str]]] = {
    "Myeloid": {
        "Macrophage": ["CD163", "MRC1", "MARCO", "CD5L", "VSIG4", "TREM2"],
        "Monocyte": ["FCN1", "CD14", "MNDA"],
    },
    "Dendritic_cell": {
        "cDC": ["CD1A", "CD1C", "CD1E", "CLEC10A", "FCER1A"],
        "pDC": ["LILRA4", "IRF8", "SPIB"],
        "Mature_DC": ["LAMP3"],
    },
    "Fibroblast": {
        "Activated_CAF": ["ACTA2", "PDGFRB", "TNC"],
        "Resting_Fibroblast": ["PDGFRA", "FBLN1", "FBN1", "ASPN", "SFRP2", "SFRP4", "DPT", "OGN"],
    },
    "Endothelial": {
        "Blood_endothelial": [
            "PECAM1",
            "VWF",
            "EGFL7",
            "CLEC14A",
            "RAMP2",
            "SOX17",
            "SOX18",
            "ERG",
        ],
        "Lymphatic_endothelial": ["PROX1", "LYVE1", "MMRN1"],
    },
    "Perivascular_SmoothMuscle": {
        "Pericyte": ["HIGD1B", "PDGFRB", "RERGL"],
        "Smooth_muscle": ["MYH11", "MYLK", "DES", "CNN1"],
    },
}
MIN_MARKERS_PER_SUBSTATE = 1


def validate_substate_markers(available_genes: set[str]) -> dict[str, dict[str, list[str]]]:
    validated: dict[str, dict[str, list[str]]] = {}
    for lineage, substates in SUBSTATE_MARKERS.items():
        validated[lineage] = {}
        for substate, genes in substates.items():
            present = [g for g in genes if g in available_genes]
            if len(present) < MIN_MARKERS_PER_SUBSTATE:
                raise PipelineError(
                    f"Substate '{lineage}/{substate}' has 0 markers present in the panel."
                )
            validated[lineage][substate] = present
    return validated


def compute_substate_scores(
    adata: ad.AnnData, layer: str, gene_pool: list[str], rng_seed: int = get_annotation_seed()
) -> pd.DataFrame:
    validated = validate_substate_markers(set(adata.var_names))
    scores = pd.DataFrame(index=adata.obs_names)
    for lineage, substates in validated.items():
        for substate, genes in substates.items():
            score_name = f"{lineage}__{substate}_score"
            sc.tl.score_genes(
                adata,
                gene_list=genes,
                gene_pool=gene_pool,
                layer=layer,
                score_name=score_name,
                random_state=rng_seed,
            )
            scores[score_name] = adata.obs[score_name]
    return scores


def build_tme_substate_inputs(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    lineage_scores_path = project_root / "data" / "derived" / "lineage_scores.parquet"
    output_path = project_root / "data" / "derived" / "tme_substate_inputs.parquet"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )
    if not lineage_scores_path.is_file():
        raise PipelineError(
            f"'{lineage_scores_path}' not found. Run `06_cell_type_annotation/02_score_major_lineages.py` first."
        )

    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    gene_pool = adata.var_names[adata.var["is_exposure_gene"]].tolist()

    scores = compute_substate_scores(adata, layer, gene_pool)

    lineage_scores = pd.read_parquet(lineage_scores_path)
    lineage_cols = [c for c in lineage_scores.columns if c.endswith("_lineage_score")]
    scores["argmax_lineage"] = (
        lineage_scores[lineage_cols].idxmax(axis=1).str.replace("_lineage_score", "", regex=False)
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(output_path)

    return {
        "n_cells": len(scores),
        "n_compartments": len(SUBSTATE_MARKERS),
        "n_cells_per_compartment": {
            lineage: int((scores["argmax_lineage"] == lineage).sum())
            for lineage in SUBSTATE_MARKERS
        },
    }
