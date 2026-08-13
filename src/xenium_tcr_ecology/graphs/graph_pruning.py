"""Tissue-gap-aware graph pruning (`09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`).

Removes Delaunay/radius edges that bridge tissue folds, holes, or gaps
between separated tissue fragments -- a well-known artefact of both graph
types (Delaunay triangulation in particular is not aware of the true
tissue footprint and will happily connect two cells across an arbitrarily
large empty region if nothing else separates them geometrically).

**Method, grounded in this project's own already-established MAD-based
outlier convention** (`MAD_FLAG_THRESHOLD = 3.5`, Iglewicz & Hoaglin 1993,
already used identically in `qc/spatial_artifacts.py`): for each section,
the Delaunay-graph edge-length distribution (`09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`'s
`delaunay.npz`) sets the reference scale for "a plausible local tissue
edge" in that section -- Delaunay, by construction, connects each cell
to its natural nearest neighbours, so its edge lengths characterise
genuine local tissue density better than an arbitrary fixed radius would.
A per-section `max_edge_length_um` is derived from this distribution's
median absolute deviation, then applied to prune both the Delaunay graph
and every radius-candidate graph for that section (a fixed-radius graph
can still bridge a gap smaller than its own radius, so pruning applies to
it too, not only to Delaunay).

**Result, checked before finalising:** on the first section
(P01_run1), the raw Delaunay edge-length distribution has median
14.6um but a dramatic outlier tail up to 1076.6um (max) -- tissue
folds/holes/fragment gaps, not measurement noise. The MAD-derived
threshold (median + 3.5*MAD/0.6745) works out to ~42.4um for that
section, consistent with the specification's tissue-scale framing
(comfortably above the ~1.5-3x median-cell-diameter scale of the smallest
candidate radii, so it does not spuriously prune genuine short-range
structure).

k-NN candidate graphs are not pruned here: the specification names
"Delaunay/radius" edges specifically, and a k-NN graph's edges are, by
construction, always to that cell's own k nearest neighbours regardless
of absolute distance -- pruning it by an absolute length threshold would
not remove "bridging" so much as disconnect isolated cells
entirely, a different problem this milestone does not attempt to solve.
Boundary-contact graphs cannot bridge a gap by construction (Phase
9.00's 0.5um tolerance already excludes anything but genuine near-contact).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

MAD_FLAG_THRESHOLD = 3.5
GRAPH_TYPES_TO_PRUNE_PREFIXES = ("delaunay", "radius_")


def modified_z_scores(values: np.ndarray) -> np.ndarray:
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros_like(values, dtype=float)
    return 0.6745 * (values - median) / mad


def compute_max_edge_length(
    edge_lengths: np.ndarray, mad_threshold: float = MAD_FLAG_THRESHOLD
) -> float:
    """Pure, testable derivation of the data-grounded pruning threshold
    from one section's Delaunay edge-length distribution."""
    if len(edge_lengths) == 0:
        return np.inf
    median = np.median(edge_lengths)
    mad = np.median(np.abs(edge_lengths - median))
    if mad == 0:
        return float(np.max(edge_lengths))
    return float(median + mad_threshold * mad / 0.6745)


def prune_long_edges(
    graph: sparse.csr_matrix, coords: np.ndarray, max_edge_length_um: float
) -> sparse.csr_matrix:
    """Pure, testable edge-length pruning for one graph."""
    coo = graph.tocoo()
    lengths = np.linalg.norm(coords[coo.row] - coords[coo.col], axis=1)
    keep = lengths <= max_edge_length_um
    pruned = sparse.csr_matrix((coo.data[keep], (coo.row[keep], coo.col[keep])), shape=graph.shape)
    return pruned


def prune_all_graphs(project_root: Path) -> dict:
    candidate_dir = project_root / "data" / "graphs" / "candidate"
    output_dir = project_root / "data" / "graphs" / "pruned"

    if not candidate_dir.is_dir():
        raise PipelineError(
            f"'{candidate_dir}' not found. Run `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py` first."
        )

    import anndata as ad

    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    adata = ad.read_h5ad(matrix_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for section_dir in sorted(d for d in candidate_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        cell_order = pd.read_csv(section_dir / "cell_order.tsv", sep="\t")["cell_id"]
        coords = adata.obs.loc[cell_order, ["x_centroid", "y_centroid"]].to_numpy()

        delaunay = sparse.load_npz(section_dir / "delaunay.npz")
        coo = delaunay.tocoo()
        upper = coo.row < coo.col
        delaunay_lengths = np.linalg.norm(coords[coo.row[upper]] - coords[coo.col[upper]], axis=1)
        max_edge_length_um = compute_max_edge_length(delaunay_lengths)

        section_out_dir = output_dir / section_id
        section_out_dir.mkdir(parents=True, exist_ok=True)
        (section_out_dir / "cell_order.tsv").write_text(
            (section_dir / "cell_order.tsv").read_text()
        )

        for graph_path in sorted(section_dir.glob("*.npz")):
            graph_name = graph_path.stem
            graph = sparse.load_npz(graph_path)
            if graph_name.startswith(GRAPH_TYPES_TO_PRUNE_PREFIXES):
                pruned = prune_long_edges(graph, coords, max_edge_length_um)
            else:
                pruned = graph
            sparse.save_npz(section_out_dir / f"{graph_name}.npz", pruned)
            manifest_rows.append(
                {
                    "section_id": section_id,
                    "graph_name": graph_name,
                    "max_edge_length_um": (
                        max_edge_length_um
                        if graph_name.startswith(GRAPH_TYPES_TO_PRUNE_PREFIXES)
                        else None
                    ),
                    "n_edges_before": int(graph.nnz / 2),
                    "n_edges_after": int(pruned.nnz / 2),
                    "n_edges_pruned": int((graph.nnz - pruned.nnz) / 2),
                }
            )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    pruned_rows = manifest[manifest["graph_name"].str.startswith(GRAPH_TYPES_TO_PRUNE_PREFIXES)]
    return {
        "n_sections": manifest["section_id"].nunique(),
        "median_max_edge_length_um": round(float(pruned_rows["max_edge_length_um"].median()), 2),
        "total_edges_pruned": int(pruned_rows["n_edges_pruned"].sum()),
        "fraction_edges_pruned_by_graph_type": (
            pruned_rows.groupby("graph_name")
            .apply(
                lambda g: round(float(g["n_edges_pruned"].sum() / g["n_edges_before"].sum()), 4),
                include_groups=False,
            )
            .to_dict()
        ),
        "output_dir": str(output_dir),
        "manifest_path": str(manifest_path),
    }
