"""Candidate spatial graph construction (`09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`).

Builds four graph types per section -- radius, k-nearest-neighbour,
Delaunay, and boundary-contact -- across a small set of candidate
parameter values for the two parametric types (radius, k-NN), not yet a
single calibrated choice: this milestone builds candidates; `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` calibrates/selects; `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`
constructs the final "primary" graph. Building several plausible
parameter values now, rather than guessing one upfront, is what makes
9.02's later calibration step meaningful.

**Candidate parameters, grounded directly in this project's own
already-established cell-geometry statistics, not arbitrary round
numbers:** `RADIUS_CANDIDATES_UM = [15.0, 30.0, 50.0]` (~1.5x, 3x, 5x the
median cell diameter of 9.82um, first established in `04_quality_control/04_estimate_transcript_spillover.py` --
spanning immediate-contact through broader-niche interaction scales).
`KNN_CANDIDATES = [6, 10, 15]` (6 and 10 already used and validated
elsewhere in this project for spatial-neighbour computations, Phase
6.06/7.02/8.04; 15 added as a broader comparison point).

**Boundary-contact graph:** cell-boundary polygons (not centroid
proxies), reusing the same `shapely.STRtree` + local/global cell-ID
reconciliation pattern already established in `04_quality_control/04_estimate_transcript_spillover.py`
. Two cells are connected
if their segmented polygons are within `BOUNDARY_CONTACT_TOLERANCE_UM`
of each other -- not a strict `predicate="intersects"` (checked
on data before finalising, not assumed: a first version using exact
intersection gave an implausible mean degree of 0.18 on a 25,964-cell
section, versus ~6 for the Delaunay graph on the same cells. Investigated
further: sampling 2,000 cell-boundary pairs, only 5.1% of each cell's
nearest neighbour is at zero distance -- Xenium's segmentation
leaves a small gap between most adjacent cells' polygons, it does not
tile the tissue seamlessly. The gap distribution has a sharp elbow:
80.5% of close pairs are captured by a 0.5um tolerance, rising only to
82.1% at 1.0um and 91.9% at 3.0um -- i.e. the genuine "touching"
population clusters tightly just above zero, and a generous tolerance
mostly just starts admitting cells that are not really adjacent).
`BOUNDARY_CONTACT_TOLERANCE_UM = 0.5` was chosen from this elbow, not
as an arbitrary round number.

Each candidate graph is saved as a symmetric scipy sparse CSR matrix
(`.npz`), aligned to that section's cell order in
`primary_analysis_matrix.h5ad` (not the raw per-section local order), so
graphs from every section share one consistent global cell-ID space and
can be indexed the same way as every other per-cell table in this
project.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import shapely
import spatialdata as sd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

RADIUS_CANDIDATES_UM = [15.0, 30.0, 50.0]
KNN_CANDIDATES = [6, 10, 15]
BOUNDARY_CONTACT_TOLERANCE_UM = 0.5


def build_radius_graph(coords: np.ndarray, radius_um: float) -> sparse.csr_matrix:
    import anndata as ad
    import squidpy as sq

    sub = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    sub.obsm["spatial"] = coords
    sq.gr.spatial_neighbors_radius(sub, radius=radius_um)
    return sub.obsp["spatial_connectivities"].tocsr()


def build_knn_graph(coords: np.ndarray, n_neighs: int) -> sparse.csr_matrix:
    import anndata as ad
    import squidpy as sq

    sub = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    sub.obsm["spatial"] = coords
    n_neighs_eff = min(n_neighs, len(coords) - 1) if len(coords) > 1 else 0
    if n_neighs_eff < 1:
        return sparse.csr_matrix((len(coords), len(coords)))
    sq.gr.spatial_neighbors_knn(sub, n_neighs=n_neighs_eff)
    return sub.obsp["spatial_connectivities"].tocsr()


def build_delaunay_graph(coords: np.ndarray) -> sparse.csr_matrix:
    import anndata as ad
    import squidpy as sq

    n = len(coords)
    if n < 4:
        return sparse.csr_matrix((n, n))
    sub = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32))
    sub.obsm["spatial"] = coords
    sq.gr.spatial_neighbors_delaunay(sub)
    return sub.obsp["spatial_connectivities"].tocsr()


def build_boundary_contact_graph(
    polygons: np.ndarray, tolerance_um: float = BOUNDARY_CONTACT_TOLERANCE_UM
) -> sparse.csr_matrix:
    """Pure, testable boundary-contact graph from an array of shapely
    polygons (already ordered to match the desired output cell order).
    Uses a small data-grounded distance tolerance
    (`predicate="dwithin"`), not strict `intersects` -- see module
    docstring for why exact intersection is far too restrictive for
    Xenium segmentation output."""
    n = len(polygons)
    if n == 0:
        return sparse.csr_matrix((0, 0))
    tree = shapely.STRtree(polygons)
    query_idx, tree_idx = tree.query(polygons, predicate="dwithin", distance=tolerance_um)
    keep = query_idx != tree_idx
    query_idx, tree_idx = query_idx[keep], tree_idx[keep]
    data = np.ones(len(query_idx), dtype=np.float32)
    return sparse.csr_matrix((data, (query_idx, tree_idx)), shape=(n, n))


def build_candidate_graphs_for_section(
    zarr_path: Path, section_cell_ids: pd.Index, coords: np.ndarray
) -> dict[str, sparse.csr_matrix]:
    """Builds every candidate graph for one section. `section_cell_ids`
    (global IDs, e.g. 'P01_run1_aaadggoi-1') and `coords` (aligned, same
    order) define the fixed output cell order every graph shares."""
    graphs: dict[str, sparse.csr_matrix] = {}
    for radius in RADIUS_CANDIDATES_UM:
        graphs[f"radius_{radius}um"] = build_radius_graph(coords, radius)
    for k in KNN_CANDIDATES:
        graphs[f"knn_{k}"] = build_knn_graph(coords, k)
    graphs["delaunay"] = build_delaunay_graph(coords)

    section_id = zarr_path.stem
    sdata = sd.read_zarr(zarr_path)
    cell_shapes = sdata["cell_boundaries"]
    prefix = f"{section_id}_"
    local_ids = [cid[len(prefix) :] if cid.startswith(prefix) else None for cid in section_cell_ids]
    available_mask = np.array([lid is not None and lid in cell_shapes.index for lid in local_ids])
    ordered_polygons = np.empty(len(section_cell_ids), dtype=object)
    if available_mask.any():
        available_local = [local_ids[i] for i in range(len(local_ids)) if available_mask[i]]
        polys = cell_shapes.loc[available_local, "geometry"].to_numpy()
        ordered_polygons[available_mask] = polys
    # Missing-boundary cells (should be rare -- `04_quality_control/03_assess_segmentation_quality.py` already
    # characterised this) get a degenerate empty point geometry so the
    # STRtree query still runs over a complete, fixed-length array;
    # an empty geometry cannot intersect anything, so such a cell
    # simply gets no boundary-contact edges, not a crash.
    for i in np.where(~available_mask)[0]:
        ordered_polygons[i] = shapely.Point()
    graphs["boundary_contact"] = build_boundary_contact_graph(ordered_polygons)

    return graphs


def build_all_candidate_graphs(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_dir = project_root / "data" / "graphs" / "candidate"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    obs = adata.obs[["section_id", "x_centroid", "y_centroid"]]

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for section_id, group in obs.groupby("section_id", observed=True):
        zarr_path = spatialdata_root / f"{section_id}.zarr"
        if not zarr_path.exists():
            raise PipelineError(
                f"'{zarr_path}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
            )

        coords = group[["x_centroid", "y_centroid"]].to_numpy()
        graphs = build_candidate_graphs_for_section(zarr_path, group.index, coords)

        section_dir = output_dir / section_id
        section_dir.mkdir(parents=True, exist_ok=True)
        cell_order_path = section_dir / "cell_order.tsv"
        pd.Series(group.index, name="cell_id").to_csv(cell_order_path, sep="\t", index=False)

        for graph_name, matrix in graphs.items():
            sparse.save_npz(section_dir / f"{graph_name}.npz", matrix)
            manifest_rows.append(
                {
                    "section_id": section_id,
                    "graph_name": graph_name,
                    "n_cells": matrix.shape[0],
                    "n_edges": int(matrix.nnz / 2),
                    "mean_degree": (
                        round(float(matrix.nnz / matrix.shape[0]), 4) if matrix.shape[0] else 0.0
                    ),
                }
            )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_sections": obs["section_id"].nunique(),
        "n_graph_types_per_section": manifest["graph_name"].nunique(),
        "n_graphs_total": len(manifest),
        "median_mean_degree_by_graph_type": manifest.groupby("graph_name")["mean_degree"]
        .median()
        .round(2)
        .to_dict(),
        "output_dir": str(output_dir),
        "manifest_path": str(manifest_path),
    }
