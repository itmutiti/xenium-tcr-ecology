"""Sender-receiver pair definitions (`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`).

Predeclares 4 biologically motivated sender-receiver comparisons
(tumour, fibroblast, myeloid, APC, each paired with T cell as receiver),
each grounded in this project's already-established findings rather
than a generic, unmotivated list:

1. `Epithelial_Tumour -> T_cell`: direct tumour-immune checkpoint
   signalling. `10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py` already found strong local depletion
   between these two lineages (z=-20.80); `11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py` already quantified
   T-cell-clone tumour engagement directly. This pair tests whether
   checkpoint/exhaustion ligand-receptor signalling specifically
   explains that spatial pattern.
2. `Fibroblast -> T_cell`: barrier signalling. `11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py` already
   found a non-trivial fibroblast barrier fraction (mean 0.201 among
   evaluable clone-sections) along the shortest path from T-cell clones
   to tumour.
3. `Myeloid -> T_cell`: activation/exhaustion-modulating signalling.
   `10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py` found `Myeloid`/`T_cell` neighbourhood enrichment
   (z=11.84); `11_clone_spatial_descriptors/03_quantify_clone_apc_support.py` found macrophage engagement enrichment
   (1.441x).
4. `Dendritic_cell -> T_cell`: antigen-presentation signalling.
   `11_clone_spatial_descriptors/03_quantify_clone_apc_support.py` found a modest positive antigen-presentation score
   excess (+0.0608) in T-cell-clone spatial neighbourhoods.

Each pair is associated with the most biologically relevant programme(s)
already validated in `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s `PROGRAM_GENE_SETS`, plus
one new `chemokine` gene set defined here (`CCL5`, `CCL19`, `CCR7`,
`CXCR4`, confirmed present in the analysis matrix against `metadata/
feature_annotation.tsv` before use; distinct from `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s `interferon` set's
`CXCL9`/`CXCL10`, which are simultaneously chemokines and
interferon-stimulated genes -- a biological overlap, not a duplication
error).

**Panel gap, flagged here at the first point in Stage 14 where it
matters:** "TGF-beta" (explicitly named in `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s scaffold, "checkpoint, chemokine,
TGF-beta, interferon and antigen-presentation programs") has no
coverage in this project's 623-gene panel -- checked against 10
canonical TGF-beta pathway genes (`TGFB1/2/3`, `TGFBR1/2/3`,
`SMAD2/3/4/7`); none are present. This programme cannot be computed
from this panel at all. `TGFB_GENE_SET` is deliberately an empty list,
not omitted, so every downstream Spatial Interactions and Barriers script can
check `len(TGFB_GENE_SET) == 0` and skip it explicitly rather than each
independently rediscovering this gap.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.program_scores import PROGRAM_GENE_SETS

CHEMOKINE_GENE_SET = ["CCL5", "CCL19", "CCR7", "CXCR4"]
TGFB_GENE_SET: list[str] = []  # checked, zero panel coverage -- see module docstring

SENDER_RECEIVER_PAIRS = [
    {
        "pair_id": "tumour_to_t_cell",
        "sender": "Epithelial_Tumour",
        "receiver": "T_cell",
        "rationale": (
            "Local depletion already found (`10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py`, z=-20.80) and clone-tumour "
            "engagement already quantified (`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`); tests whether checkpoint signalling "
            "explains the spatial exclusion pattern."
        ),
        "relevant_programs": ["exhaustion", "antigen_presentation"],
    },
    {
        "pair_id": "fibroblast_to_t_cell",
        "sender": "Fibroblast",
        "receiver": "T_cell",
        "rationale": (
            "Non-trivial fibroblast barrier fraction already found along the shortest "
            "path from T-cell clones to tumour (`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`, mean 0.201 among evaluable "
            "clone-sections)."
        ),
        "relevant_programs": ["chemokine", "tgf_beta"],
    },
    {
        "pair_id": "myeloid_to_t_cell",
        "sender": "Myeloid",
        "receiver": "T_cell",
        "rationale": (
            "Myeloid/T_cell neighbourhood enrichment already found (`10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py`, z=11.84) "
            "and macrophage engagement enrichment already found (`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`, 1.441x)."
        ),
        "relevant_programs": ["activation", "interferon", "chemokine"],
    },
    {
        "pair_id": "dendritic_cell_to_t_cell",
        "sender": "Dendritic_cell",
        "receiver": "T_cell",
        "rationale": (
            "Modest positive antigen-presentation score excess already found in T-cell-clone "
            "spatial neighbourhoods (`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`, +0.0608)."
        ),
        "relevant_programs": ["antigen_presentation", "activation"],
    },
]


def build_program_gene_sets_for_interactions() -> dict[str, list[str]]:
    """`05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s validated programme gene sets
    extended with Spatial Interactions and Barriers's new `chemokine` set and the
    empty `tgf_beta` set (zero panel coverage, not omitted)."""
    extended = dict(PROGRAM_GENE_SETS)
    extended["chemokine"] = CHEMOKINE_GENE_SET
    extended["tgf_beta"] = TGFB_GENE_SET
    return extended


def validate_sender_receiver_pairs(
    pairs: list[dict], valid_lineages: set[str], valid_programs: set[str]
) -> list[str]:
    """Validation errors for a candidate pair list -- unknown
    sender/receiver lineage or unknown programme name. Returns an empty
    list if every pair is valid."""
    errors = []
    for pair in pairs:
        if pair["sender"] not in valid_lineages:
            errors.append(f"{pair['pair_id']}: unknown sender lineage '{pair['sender']}'")
        if pair["receiver"] not in valid_lineages:
            errors.append(f"{pair['pair_id']}: unknown receiver lineage '{pair['receiver']}'")
        for program in pair["relevant_programs"]:
            if program not in valid_programs:
                errors.append(f"{pair['pair_id']}: unknown programme '{program}'")
    return errors


def build_sender_receiver_config(project_root: Path) -> dict:
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    output_path = project_root / "config" / "sender_receiver_pairs.yaml"

    if not final_annotations_path.is_file():
        raise PipelineError(
            f"'{final_annotations_path}' not found. Run `06_cell_type_annotation/06_integrate_annotation_evidence.py` first."
        )

    import pandas as pd

    final_annotations = pd.read_parquet(final_annotations_path)
    valid_lineages = set(final_annotations["final_lineage"].dropna().unique())

    program_gene_sets = build_program_gene_sets_for_interactions()
    valid_programs = set(program_gene_sets.keys())

    errors = validate_sender_receiver_pairs(SENDER_RECEIVER_PAIRS, valid_lineages, valid_programs)
    if errors:
        raise PipelineError(f"Validation error(s) in SENDER_RECEIVER_PAIRS: {errors}")

    config = {
        "sender_receiver_pairs": SENDER_RECEIVER_PAIRS,
        "program_gene_sets": program_gene_sets,
        "programs_with_zero_panel_coverage": [
            name for name, genes in program_gene_sets.items() if len(genes) == 0
        ],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(config, sort_keys=False, default_flow_style=False))

    return {
        "n_pairs": len(SENDER_RECEIVER_PAIRS),
        "n_programs": len(program_gene_sets),
        "programs_with_zero_coverage": config["programs_with_zero_panel_coverage"],
        "output_path": str(output_path),
    }
