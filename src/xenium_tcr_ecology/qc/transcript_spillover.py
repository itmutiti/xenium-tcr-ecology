"""Transcript spillover risk estimation (`04_quality_control/04_estimate_transcript_spillover.py`).

Xenium (and imaging-based spatial transcriptomics generally) is subject to a
known, near-field artefact: transcripts emitted near a cell's segmentation
boundary can be assigned to the wrong cell when two cells sit close
together, especially across a predicted-cell-type boundary (e.g. an
immune cell pressed against a tumour cell) where a spillover transcript
would masquerade as evidence of the wrong cell's identity/state. This phase
does not correct counts -- it flags, per cell, how much a cell's expression
profile should be trusted given its physical neighbourhood, using the
predicted cell type now available from Cell Type Annotation, since this
milestone was deferred until after Cell Type Annotation.

Design, grounded in cell geometry (not centroid-distance proxies) and a
data-derived search radius:

  - SEARCH_RADIUS_UM = 10.0 -- approximately one median cell diameter.
    Checked against `cell_qc_metrics.parquet`'s `cell_area` column
    for all 1,186,916 cells: `2*sqrt(area/pi)` has median 9.82um (IQR
    8.05-11.89um). Ten microns is therefore a defensible "immediately
    adjacent cell" scale -- transcript spillover is a near-field effect
    between physically close/touching cells, not a tissue-region-scale
    effect, so a substantially larger radius would dilute the signal with
    cells that are not plausible spillover sources.
  - For each cell, cell-boundary polygons (not centroids) are used via
    a vectorised `shapely.STRtree` `dwithin` query, matching this project's
    established preference (`04_quality_control/03_assess_segmentation_quality.py`) for exact polygon geometry over
    approximations. Edge-to-edge polygon distance is exact: touching or
    overlapping polygons have distance 0.
  - `spillover_risk_score`: for each cell, the proximity-weighted fraction
    of its physical neighbours (within SEARCH_RADIUS_UM) that carry a
    different predicted major lineage (`06_cell_type_annotation/06_integrate_annotation_evidence.py`'s `final_lineage`) --
    `mean over neighbours of (1 - distance/SEARCH_RADIUS_UM) *
    1[different lineage]`. This directly implements the specification's
    "distance-to-boundary and neighbour-identity weighting": closer
    different-type neighbours contribute more; same-type neighbours
    contribute nothing; a cell with no nearby neighbours at all scores 0
    (nothing to spill from). Bounded in [0, 1] by construction.
  - `is_boundary_adjacent_to_different_type`: a stricter boolean flag for
    the literal case the specification describes -- at least one
    touching/overlapping (distance <= TOUCHING_EPSILON_UM)
    neighbour of a different predicted lineage. This is the highest-
    confidence spillover-risk indicator, independent of the continuous
    score's proximity weighting.

Scope decision: both the target cell and candidate neighbours are drawn
from `final_cell_annotations.parquet` (`06_cell_type_annotation/06_integrate_annotation_evidence.py`'s QC-passed, typed
cohort), not the raw pre-QC cell set -- a neighbour's identity is only
informative if it has a predicted type, and QC-excluded cells already carry
no trustworthy expression profile of their own to have "spilled" in a
meaningful sense. This does mean a cell that is physically adjacent only to
QC-excluded neighbours is scored using its remaining QC-passed neighbours
only (or 0 if none), which is stated here rather than left implicit.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import shapely
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Grounded directly in data (see module docstring): median cell diameter
# across all 1,186,916 cells is 9.82um.
SEARCH_RADIUS_UM = 10.0

# Polygon-to-polygon distance is a continuous float; two
# touching/overlapping Xenium polygons measure exactly 0.0 in practice, but
# a small epsilon avoids any floating-point edge case being missed.
TOUCHING_EPSILON_UM = 1e-6


def compute_section_spillover_risk(
    cell_polygons: pd.Series,
    lineages: pd.Series,
    search_radius_um: float = SEARCH_RADIUS_UM,
) -> pd.DataFrame:
    """Pure geometric computation for one section's worth of cells --
    factored out so the core logic is testable with plain shapely polygons,
    not only via a full SpatialData store. `cell_polygons` and `lineages`
    must share the same index (cell IDs) and cover exactly the same cells.
    """
    if not cell_polygons.index.equals(lineages.index):
        raise PipelineError("cell_polygons and lineages must share the same index.")

    n = len(cell_polygons)
    if n == 0:
        return pd.DataFrame(
            columns=[
                "n_neighbors_within_radius",
                "n_different_type_neighbors_within_radius",
                "nearest_different_type_distance_um",
                "spillover_risk_score",
                "is_boundary_adjacent_to_different_type",
            ]
        )

    cell_ids = cell_polygons.index.to_numpy()
    polys = cell_polygons.to_numpy()
    lineage_arr = lineages.to_numpy()

    tree = shapely.STRtree(polys)
    query_idx, tree_idx = tree.query(polys, predicate="dwithin", distance=search_radius_um)

    # Drop trivial self-pairs (every polygon is dwithin itself at distance 0).
    keep = query_idx != tree_idx
    query_idx, tree_idx = query_idx[keep], tree_idx[keep]

    if len(query_idx) == 0:
        distances = np.array([])
    else:
        distances = shapely.distance(polys[query_idx], polys[tree_idx])

    is_different_type = lineage_arr[query_idx] != lineage_arr[tree_idx]
    proximity_weight = 1.0 - (distances / search_radius_um)

    pairs = pd.DataFrame(
        {
            "query_idx": query_idx,
            "distance": distances,
            "is_different_type": is_different_type,
            "weighted_contribution": np.where(is_different_type, proximity_weight, 0.0),
        }
    )

    n_neighbors = pairs.groupby("query_idx").size()
    n_different = pairs.groupby("query_idx")["is_different_type"].sum()
    weighted_sum = pairs.groupby("query_idx")["weighted_contribution"].sum()
    different_pairs = pairs[pairs["is_different_type"]]
    nearest_different_distance = different_pairs.groupby("query_idx")["distance"].min()
    is_touching_different = (
        different_pairs[different_pairs["distance"] <= TOUCHING_EPSILON_UM]
        .groupby("query_idx")
        .size()
    )

    result = pd.DataFrame(index=pd.RangeIndex(n))
    result["n_neighbors_within_radius"] = n_neighbors.reindex(result.index, fill_value=0).astype(
        int
    )
    result["n_different_type_neighbors_within_radius"] = n_different.reindex(
        result.index, fill_value=0
    ).astype(int)
    result["nearest_different_type_distance_um"] = nearest_different_distance.reindex(result.index)
    result["spillover_risk_score"] = (
        weighted_sum.reindex(result.index, fill_value=0.0)
        / result["n_neighbors_within_radius"].replace(0, np.nan)
    ).fillna(0.0)
    result["is_boundary_adjacent_to_different_type"] = (
        is_touching_different.reindex(result.index, fill_value=0) > 0
    )
    result.index = cell_ids
    result.index.name = "cell_id"
    return result


def build_spillover_risk_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_path = project_root / "data" / "derived" / "spillover_risk.parquet"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )
    if not final_annotations_path.is_file():
        raise PipelineError(
            f"'{final_annotations_path}' not found. Run `06_cell_type_annotation/06_integrate_annotation_evidence.py` first."
        )

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    final_annotations = pd.read_parquet(final_annotations_path)
    section_ids = adata.obs["section_id"].reindex(final_annotations.index)

    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        raise PipelineError(
            f"No .zarr stores found under '{spatialdata_root}'. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    section_results = []
    n_cells_missing_polygon = 0
    for zarr_path in zarr_paths:
        section_id = zarr_path.stem
        section_cell_ids = section_ids.index[section_ids == section_id]
        if len(section_cell_ids) == 0:
            continue

        sdata = sd.read_zarr(zarr_path)
        cell_shapes = sdata["cell_boundaries"]

        # `03_spatialdata_import/05_build_combined_analysis_object.py`'s combine_sections.py globally-uniquifies obs_names as
        # f"{section_id}_{local_cell_id}" (necessary since Xenium's raw
        # per-section cell IDs collide across sections). Per-section zarr
        # stores retain the original local (unprefixed) IDs, so the two ID
        # spaces must be reconciled here before any join.
        prefix = f"{section_id}_"
        local_to_global = {
            cid[len(prefix) :]: cid for cid in section_cell_ids if cid.startswith(prefix)
        }
        available_local = cell_shapes.index.intersection(local_to_global.keys())
        n_cells_missing_polygon += len(section_cell_ids) - len(available_local)
        if len(available_local) == 0:
            continue

        available_global = [local_to_global[lid] for lid in available_local]
        cell_polygons = cell_shapes.loc[available_local, "geometry"]
        cell_polygons.index = available_global
        lineages = final_annotations.loc[available_global, "final_lineage"]
        section_result = compute_section_spillover_risk(cell_polygons, lineages)
        section_results.append(section_result)

    if not section_results:
        raise PipelineError("No section produced any spillover-risk results.")

    result = pd.concat(section_results)
    result = result.join(final_annotations[["final_lineage"]])
    result["section_id"] = section_ids.reindex(result.index)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    risk_by_lineage = (
        result.groupby("final_lineage", observed=True)["spillover_risk_score"]
        .mean()
        .round(4)
        .to_dict()
    )
    return {
        "n_cells": len(result),
        "n_sections": result["section_id"].nunique(),
        "n_cells_missing_polygon": n_cells_missing_polygon,
        "mean_spillover_risk_score": round(float(result["spillover_risk_score"].mean()), 4),
        "fraction_boundary_adjacent_to_different_type": round(
            float(result["is_boundary_adjacent_to_different_type"].mean()), 4
        ),
        "mean_spillover_risk_score_by_lineage": risk_by_lineage,
        "search_radius_um": SEARCH_RADIUS_UM,
    }
