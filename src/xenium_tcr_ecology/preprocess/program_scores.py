"""Curated program score calculation (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`).

Computes per-cell scores for 8 named biological programs using scanpy's
standard `score_genes` method (mean expression of program genes minus a
random, expression-matched control gene set -- the same approach
underlying Seurat's `AddModuleScore`), on the primary normalisation layer
selected in `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` (`adata.uns["primary_normalization_layer"]`, read
programmatically rather than hardcoded, so this script always uses
whichever layer `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` actually decided on).

Every marker gene below was checked against the 399-gene
`biological_gene` panel (`metadata/feature_annotation.tsv`), not assumed
from a generic literature list -- this is a targeted, immuno-oncology/TCR
-focused panel (heavy T-cell/NK/myeloid coverage, consistent with the
earlier marker-coverage check in Quality Control for a different
purpose), and several programs have load-bearing gene-set gaps relative
to their canonical literature signatures:

  - cytotoxicity, exhaustion, proliferation: strong coverage, close to
    canonical signatures.
  - activation: reasonable coverage (5 genes), missing some canonical
    markers (ICOS, TNFRSF4/OX40, CD40LG) not in this panel.
  - interferon: thin (3 genes). No canonical interferon-stimulated genes
    (STAT1, IRF1, ISG15, MX1, OAS1, IFIT1/3) are in this panel at all;
    CXCL9/CXCL10 (IFN-gamma-induced chemokines) and IRF8 are used as a
    proxy, not a validated ISG signature.
  - stress: thin and non-canonical. No classical heat-shock/immediate-
    early genes (HSPA1A/B, FOS, JUN, DNAJB1, HSPB1) are in this panel;
    GDF15/NFE2L2/KEAP1/STC1/STC2 (oxidative/metabolic stress response) are
    used instead -- a different biological axis from the
    "cellular stress response" this program name usually implies.
  - EMT: thin and non-canonical. No canonical EMT transcription factors
    (SNAI2, ZEB1, TWIST1, CDH2) or VIM/FN1 are in this panel; the set used
    here (SNAI1 plus stromal/mesenchymal markers ACTA2, PDGFRA, PDGFRB,
    THY1, TNC, THBS2, VCAN) is closer to a general
    "mesenchymal/stromal-associated" signature than a validated
    epithelial-to-mesenchymal *transition* signature.
  - antigen_presentation: thin and MHC-II/DC-skewed. No MHC class I genes
    (HLA-A/B/C, B2M, TAP1/2) are in this panel at all -- a materially
    important gap given tumour immune-evasion narratives are often about
    MHC-I loss specifically. The set used here (HLA-DRA, HLA-DQB2, CD1A/C/E,
    LAMP3) captures MHC-II / dendritic-cell antigen presentation only.

These gaps are inherent to the panel design -- this panel prioritises
T-cell/TCR coverage over broad tumour/stress/EMT biology -- and are recorded
here, not silently absorbed into a generically-named score.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

PROGRAM_GENE_SETS: dict[str, list[str]] = {
    "cytotoxicity": [
        "GZMA",
        "GZMB",
        "GZMK",
        "PRF1",
        "GNLY",
        "NKG7",
        "KLRD1",
        "KLRB1",
        "KLRC1",
        "FGFBP2",
    ],
    "exhaustion": ["PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1", "TOX"],
    "activation": ["CD69", "IL2RA", "TNFRSF9", "CD28", "CD27"],
    "interferon": ["CXCL9", "CXCL10", "IRF8"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "CCNB2", "CDK1", "UBE2C", "CENPF"],
    "stress": ["GDF15", "NFE2L2", "KEAP1", "STC1", "STC2"],
    "emt": ["SNAI1", "ACTA2", "PDGFRA", "PDGFRB", "THY1", "TNC", "THBS2", "VCAN"],
    "antigen_presentation": ["HLA-DRA", "HLA-DQB2", "CD1A", "CD1C", "CD1E", "LAMP3"],
}

# Programs whose gene set diverges meaningfully from the canonical
# literature signature due to panel gaps (see module docstring) --
# surfaced programmatically so downstream reports/scripts don't need to
# re-derive this from prose.
THIN_COVERAGE_PROGRAMS = {"interferon", "stress", "emt", "antigen_presentation"}

MIN_GENES_PER_PROGRAM = 2


def validate_gene_sets(
    available_genes: set[str], gene_sets: dict[str, list[str]] = PROGRAM_GENE_SETS
) -> dict[str, list[str]]:
    """Restricts each program's gene set to genes actually present in
    `available_genes`, raising if any program would be left with too few
    genes to compute a meaningful score. Pure and testable independently of
    scanpy/AnnData."""
    validated = {}
    for program, genes in gene_sets.items():
        present = [g for g in genes if g in available_genes]
        if len(present) < MIN_GENES_PER_PROGRAM:
            raise PipelineError(
                f"Program '{program}' has only {len(present)} gene(s) present in the panel "
                f"(of {len(genes)} curated) -- below the minimum of {MIN_GENES_PER_PROGRAM} "
                "needed for a meaningful score."
            )
        validated[program] = present
    return validated


def compute_program_scores(
    adata: ad.AnnData, layer: str, gene_pool: list[str], rng_seed: int = get_annotation_seed()
) -> pd.DataFrame:
    validated_sets = validate_gene_sets(set(adata.var_names))

    scores = pd.DataFrame(index=adata.obs_names)
    for program, genes in validated_sets.items():
        score_name = f"{program}_score"
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


def build_program_scores_report(project_root: Path) -> dict:
    analysis_ready_path = project_root / "data" / "objects" / "analysis_ready.h5ad"
    output_path = project_root / "data" / "derived" / "program_scores.parquet"

    if not analysis_ready_path.is_file():
        raise PipelineError(
            f"'{analysis_ready_path}' not found. Run `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py` first."
        )

    adata = ad.read_h5ad(analysis_ready_path)
    if "primary_normalization_layer" not in adata.uns:
        raise PipelineError(
            f"'{analysis_ready_path}' has no uns['primary_normalization_layer'] -- run `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` first "
            "."
        )
    layer = adata.uns["primary_normalization_layer"]
    if layer not in adata.layers:
        raise PipelineError(
            f"Primary normalisation layer '{layer}' not found in '{analysis_ready_path}'.layers."
        )

    gene_pool = adata.var_names[adata.var["is_exposure_gene"]].tolist()

    scores = compute_program_scores(adata, layer=layer, gene_pool=gene_pool)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(output_path)

    validated_sets = validate_gene_sets(set(adata.var_names))
    return {
        "n_cells": len(scores),
        "n_programs": len(validated_sets),
        "layer_used": layer,
        "genes_per_program": {k: len(v) for k, v in validated_sets.items()},
        "thin_coverage_programs": sorted(THIN_COVERAGE_PROGRAMS & set(validated_sets.keys())),
    }
