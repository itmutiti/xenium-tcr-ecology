"""Clone antigen-presenting-cell (APC) support (`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`).

Three proximity/enrichment quantities per (clone, section) -- Phase
11.00's established unit -- each reusing existing infrastructure rather
than a fresh, independently-invented computation:

1. **Dendritic cell adjacency**: reused directly from `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`'s
   `radius_30.0um__Dendritic_cell` composition column, exactly as Phase
   11.02 reused the `Epithelial_Tumour` column for malignant adjacency.
2. **Macrophage adjacency**: `Macrophage` is a `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R` myeloid
   substate, not one of `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`'s 12 major-lineage composition
   columns, so no existing column covers it -- computed here via Phase
   10.01's `compute_composition_vectors` function (`xenium_tcr_
   ecology.niches.local_composition`), reused directly (not
   reimplemented) with a binary Macrophage/non-Macrophage indicator in
   place of the full 12-lineage one-hot, on the same calibrated Phase
   9.03 30um primary graph.
3. **Antigen-presentation programme proximity**: Preprocessing and Normalisation's already-scored
   `antigen_presentation_score` (`data/derived/program_scores.
   parquet`), aggregated over each cell's spatial neighbours via the
   same graph-multiplication technique as (1)/(2) but for a continuous
   score rather than a categorical composition
   (`compute_neighbourhood_mean_score`, a small, analogous function) --
   a clone cell's own score is not what matters here (a T cell does not
   itself run an antigen-presentation programme); what matters is
   whether antigen-presentation activity is elevated in its spatial
   neighbourhood.

Dendritic-cell and macrophage adjacency are bounded [0, 1] composition
fractions, opportunity-normalised as a `*_engagement_ratio` (clone's
neighbourhood value / that section's mean), the identical convention
established in `11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`. Antigen-presentation
score is different in kind -- a centred, potentially-negative z-scored
gene-programme score (checked against `program_scores.
parquet`: mean ~-0.015, range roughly -0.44 to +2.12), for which a
ratio is mathematically nonsensical (division by a near-zero or
negative baseline can flip sign or blow up arbitrarily -- caught on the
first run, which produced an uninterpretable mean ratio of -10.27).
Normalised instead as `antigen_presentation_score_excess`, a difference
(clone value minus section baseline, `compute_score_excess`) -- the
correct opportunity-normalisation for a signed, centred score.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.clone_ecology.spatial_descriptors import filter_to_primary_cohort
from xenium_tcr_ecology.clone_ecology.tumour_engagement import compute_engagement_ratio
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.niches.local_composition import compute_composition_vectors

DC_ADJACENCY_COLUMN = "radius_30.0um__Dendritic_cell"
MACROPHAGE_SUBSTATE = "Macrophage"


def compute_score_excess(clone_value: float, section_baseline: float) -> float:
    """Pure, testable: opportunity-normalised difference (not ratio) for
    a centred, potentially-negative continuous score.
    `compute_engagement_ratio` (`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`) is only valid for bounded
    [0, 1] proportions -- `antigen_presentation_score` is a z-scored
    gene-programme score (mean ~-0.015, range roughly -0.44 to +2.12;
    checked against `program_scores.parquet` before choosing
    this), so a ratio against a near-zero or negative baseline is
    mathematically nonsensical (division can flip sign or blow up
    arbitrarily) -- a difference is the correct normalisation."""
    if section_baseline is None or np.isnan(section_baseline):
        return float("nan")
    return float(clone_value - section_baseline)


def compute_neighbourhood_mean_score(graph: sparse.csr_matrix, score: pd.Series) -> pd.Series:
    """Pure, testable: graph-neighbour mean of a continuous per-cell
    score (row-normalised by degree), NaN for zero-degree cells -- the
    same convention as `compute_composition_vectors`, generalised from a
    categorical composition to a continuous score."""
    values = score.to_numpy(dtype=np.float64)
    weighted_sum = np.asarray(graph @ values).ravel()
    degree = np.asarray(graph.sum(axis=1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_score = weighted_sum / degree
    mean_score[degree == 0] = np.nan
    return pd.Series(mean_score, index=score.index)


def build_clone_apc_support(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    local_compositions_path = project_root / "data" / "derived" / "local_compositions.parquet"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    program_scores_path = project_root / "data" / "derived" / "program_scores.parquet"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    output_path = project_root / "data" / "derived" / "clone_apc_support.parquet"

    for path, phase in [
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
        (sample_manifest_path, None),
        (
            local_compositions_path,
            "`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`",
        ),
        (final_annotations_path, "`06_cell_type_annotation/06_integrate_annotation_evidence.py`"),
        (program_scores_path, "Preprocessing and Normalisation"),
        (
            primary_graphs_dir,
            "`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found." + (f" Run {phase} first." if phase else ""))

    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
    primary_sections = set(filter_to_primary_cohort(sample_manifest)["section_id"])

    resolved = pd.read_parquet(resolved_calls_path)
    clonal = resolved[resolved["resolution"].isin(["singlet", "low_confidence"])].copy()
    clonal["clone_id"] = clonal["detected_probes"]

    high_confidence_clones = pd.read_parquet(high_confidence_clones_path)
    release_clone_ids = set(high_confidence_clones["clone_id"])
    clonal = clonal[
        clonal["clone_id"].isin(release_clone_ids) & clonal["section_id"].isin(primary_sections)
    ]
    if len(clonal) == 0:
        raise PipelineError(
            "No primary-cohort, release clone cells found -- cannot build APC support."
        )

    local_compositions = pd.read_parquet(local_compositions_path)
    final_annotations = pd.read_parquet(final_annotations_path)
    program_scores = pd.read_parquet(program_scores_path)
    is_macrophage = final_annotations["final_substate"] == MACROPHAGE_SUBSTATE

    macrophage_adjacency_parts = []
    antigen_presentation_parts = []
    for section_id in sorted(set(clonal["section_id"].unique())):
        node_metadata = pd.read_csv(
            primary_graphs_dir / section_id / "node_metadata.tsv", sep="\t"
        ).set_index("cell_id")
        graph = sparse.load_npz(primary_graphs_dir / section_id / "primary_graph.npz")

        macrophage_indicator = is_macrophage.reindex(node_metadata.index).map(
            {True: MACROPHAGE_SUBSTATE, False: "Other"}
        )
        macrophage_composition = compute_composition_vectors(
            graph, macrophage_indicator, [MACROPHAGE_SUBSTATE]
        )
        macrophage_adjacency_parts.append(macrophage_composition[MACROPHAGE_SUBSTATE])

        section_scores = program_scores.reindex(node_metadata.index)["antigen_presentation_score"]
        antigen_presentation_parts.append(compute_neighbourhood_mean_score(graph, section_scores))

    macrophage_adjacency = pd.concat(macrophage_adjacency_parts)
    antigen_presentation_neighbourhood = pd.concat(antigen_presentation_parts)

    section_dc_mean = local_compositions.groupby("section_id")[DC_ADJACENCY_COLUMN].mean()
    section_macrophage_mean = macrophage_adjacency.groupby(
        local_compositions["section_id"].reindex(macrophage_adjacency.index)
    ).mean()
    section_antigen_presentation_mean = antigen_presentation_neighbourhood.groupby(
        local_compositions["section_id"].reindex(antigen_presentation_neighbourhood.index)
    ).mean()

    rows = []
    for (clone_id, section_id), group in clonal.groupby(["clone_id", "section_id"], observed=True):
        cell_ids = pd.Index(group.index)

        clone_dc = local_compositions.reindex(cell_ids)[DC_ADJACENCY_COLUMN]
        dc_adjacency = float(clone_dc.mean()) if clone_dc.notna().any() else float("nan")
        dc_engagement_ratio = compute_engagement_ratio(
            dc_adjacency, section_dc_mean.get(section_id, float("nan"))
        )

        clone_macrophage = macrophage_adjacency.reindex(cell_ids)
        macrophage_adjacency_value = (
            float(clone_macrophage.mean()) if clone_macrophage.notna().any() else float("nan")
        )
        macrophage_engagement_ratio = compute_engagement_ratio(
            macrophage_adjacency_value, section_macrophage_mean.get(section_id, float("nan"))
        )

        clone_antigen_presentation = antigen_presentation_neighbourhood.reindex(cell_ids)
        antigen_presentation_value = (
            float(clone_antigen_presentation.mean())
            if clone_antigen_presentation.notna().any()
            else float("nan")
        )
        antigen_presentation_score_excess = compute_score_excess(
            antigen_presentation_value,
            section_antigen_presentation_mean.get(section_id, float("nan")),
        )

        rows.append(
            {
                "clone_id": clone_id,
                "section_id": section_id,
                "patient_id": group["patient_id"].iloc[0],
                "n_cells": len(group),
                "dc_adjacency": dc_adjacency,
                "dc_engagement_ratio": dc_engagement_ratio,
                "macrophage_adjacency": macrophage_adjacency_value,
                "macrophage_engagement_ratio": macrophage_engagement_ratio,
                "antigen_presentation_neighbourhood_score": antigen_presentation_value,
                "antigen_presentation_score_excess": antigen_presentation_score_excess,
            }
        )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_clone_section_rows": len(result),
        "n_distinct_clones": int(result["clone_id"].nunique()),
        "mean_dc_engagement_ratio": float(result["dc_engagement_ratio"].mean()),
        "mean_macrophage_engagement_ratio": float(result["macrophage_engagement_ratio"].mean()),
        "mean_antigen_presentation_score_excess": float(
            result["antigen_presentation_score_excess"].mean()
        ),
        "output_path": str(output_path),
    }
