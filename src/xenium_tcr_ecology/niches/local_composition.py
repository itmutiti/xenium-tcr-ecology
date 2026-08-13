"""Local neighbourhood composition vectors (`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`).

For every cell, computes the fraction of each of `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s 12 major
lineages among its spatial neighbours (`09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`'s gap-pruned
graphs), at the same three candidate scales already established and
validated in `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`, `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py` (radius 15/30/50um) -- not a new,
independently-chosen scale set, so this output is directly comparable to
everything already computed at those same scales (`09_spatial_graph_construction_and_calibration/06_run_graph_sensitivity_grid.py`'s
sensitivity grid in particular).

**Method:** an efficient sparse-matrix computation, not a per-cell
Python loop -- neighbour lineage counts are `graph @ one_hot_lineage`
(one sparse matrix multiply per section per scale), then row-normalised
by each cell's degree at that scale to give fractions. A cell's own
lineage is excluded from its own composition vector (standard cellular-
neighbourhood convention, e.g. Schurch et al. 2020's CN method) -- the
composition describes what surrounds a cell, not the cell itself.

Cells with zero neighbours at a given scale (an already-characterised
possibility -- `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` found `radius_15um` leaves many cells poorly
connected) get an all-NaN composition vector at that scale, not a
fabricated all-zero one -- "no neighbours to describe" is a
different condition from "neighbours exist and are evenly mixed."
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

SCALES = ["radius_15.0um", "radius_30.0um", "radius_50.0um"]
LINEAGE_COLUMN = "final_lineage"


def compute_composition_vectors(
    graph: sparse.csr_matrix, lineage: pd.Series, lineages_order: list[str]
) -> pd.DataFrame:
    """Pure, testable neighbour-composition computation for one section
    at one scale. Returns a DataFrame indexed like `lineage`, one column
    per lineage in `lineages_order`, all-NaN rows for zero-degree cells."""
    n = len(lineage)
    one_hot = np.zeros((n, len(lineages_order)), dtype=np.float64)
    lineage_to_col = {name: i for i, name in enumerate(lineages_order)}
    for i, value in enumerate(lineage.to_numpy()):
        col = lineage_to_col.get(value)
        if col is not None:
            one_hot[i, col] = 1.0

    neighbour_counts = graph @ one_hot
    degree = np.asarray(graph.sum(axis=1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        composition = neighbour_counts / degree[:, None]
    composition[degree == 0, :] = np.nan

    return pd.DataFrame(composition, index=lineage.index, columns=lineages_order)


def build_local_compositions(project_root: Path) -> dict:
    pruned_dir = project_root / "data" / "graphs" / "pruned"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    output_path = project_root / "data" / "derived" / "local_compositions.parquet"

    if not pruned_dir.is_dir():
        raise PipelineError(
            f"'{pruned_dir}' not found. Run `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py` first."
        )
    if not final_annotations_path.is_file():
        raise PipelineError(
            f"'{final_annotations_path}' not found. Run `06_cell_type_annotation/06_integrate_annotation_evidence.py` first."
        )

    final_annotations = pd.read_parquet(final_annotations_path)
    lineages_order = sorted(final_annotations[LINEAGE_COLUMN].dropna().unique())

    section_results = []
    n_cells_zero_degree_by_scale = {scale: 0 for scale in SCALES}
    for section_dir in sorted(d for d in pruned_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        cell_order = pd.read_csv(section_dir / "cell_order.tsv", sep="\t")["cell_id"]
        lineage = final_annotations.reindex(cell_order)[LINEAGE_COLUMN]

        section_frames = []
        for scale in SCALES:
            graph_path = section_dir / f"{scale}.npz"
            if not graph_path.is_file():
                raise PipelineError(
                    f"'{graph_path}' not found. Run `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`, `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py` first."
                )
            graph = sparse.load_npz(graph_path)
            composition = compute_composition_vectors(graph, lineage, lineages_order)
            composition.columns = [f"{scale}__{c}" for c in composition.columns]
            n_cells_zero_degree_by_scale[scale] += int(composition.isna().any(axis=1).sum())
            section_frames.append(composition)

        section_result = pd.concat(section_frames, axis=1)
        section_result["section_id"] = section_id
        section_results.append(section_result)

    result = pd.concat(section_results)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_cells": len(result),
        "n_lineages": len(lineages_order),
        "n_scales": len(SCALES),
        "n_composition_columns": len(SCALES) * len(lineages_order),
        "n_cells_zero_degree_by_scale": n_cells_zero_degree_by_scale,
        "output_path": str(output_path),
    }
