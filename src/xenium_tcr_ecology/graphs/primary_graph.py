"""Primary cell graph construction (`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`).

Creates the prespecified primary graph -- one canonical graph per
section, not a further candidate -- with node and edge metadata, built
from `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`'s calibrated radius (`config/graph_parameters.yaml`,
`calibrated_radius_um`) applied to `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`'s already gap-pruned graph
for that radius.

**Judgment call, made and documented: radius, not k-NN, is the
primary graph type.** `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` calibrated both a radius
(30um) and a k (6), but the specification asks for "the" primary graph
(singular). A fixed-radius graph has a constant, biologically literal
physical-distance interpretation ("these two cells are within 30um of
each other") that does not change with local tissue density, unlike a
k-NN graph, whose edges mean something different in dense vs. sparse
regions of the same section -- and this project's later spatial-biology
questions (barrier topology, clone-tumour engagement distance) are
inherently about literal physical proximity, not a locally-adaptive
neighbour count. The calibrated k-NN graph remains available (`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`'s
own output, `data/graphs/pruned/*/knn_6.npz`) for any downstream use that
specifically wants density-adaptive neighbourhoods instead.

**Edge metadata:** the primary graph is upgraded from `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`-9.01's
binary adjacency to a weighted graph, Euclidean distance (um) as the
edge weight -- the most direct, compact way to carry edge metadata
without a separate large edge table.

**Node metadata:** joined directly from `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s
`final_cell_annotations.parquet` (`final_lineage`, `final_substate`,
`confidence`, `is_ambiguous`) plus `patient_id`/`section_id`/coordinates
from the primary analysis matrix -- every node in the primary graph
carries the same cell-type/confidence information already established
for that cell everywhere else in this project, not a graph-specific
re-derivation.

**Patient-separated components:** enforced by construction, not merely
observed -- the primary graph is built and stored per section (never
spanning sections, let alone patients), consistent with every other
per-section graph artefact in this phase. A cross-check is run and
reported (not assumed): no primary graph should have any edge connecting
two different patients, since sections themselves are already patient-
scoped 1:1.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

NODE_METADATA_COLUMNS = ["final_lineage", "final_substate", "confidence", "is_ambiguous"]


def build_weighted_graph(adjacency: sparse.csr_matrix, coords: np.ndarray) -> sparse.csr_matrix:
    """Pure, testable conversion of a binary adjacency matrix to a
    distance-weighted graph."""
    coo = adjacency.tocoo()
    distances = np.linalg.norm(coords[coo.row] - coords[coo.col], axis=1)
    return sparse.csr_matrix((distances, (coo.row, coo.col)), shape=adjacency.shape)


def build_primary_graphs(project_root: Path) -> dict:
    pruned_dir = project_root / "data" / "graphs" / "pruned"
    graph_params_path = project_root / "config" / "graph_parameters.yaml"
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    output_dir = project_root / "data" / "graphs" / "primary_graphs"

    for p in (pruned_dir, graph_params_path, matrix_path, final_annotations_path):
        if not p.is_file() and not p.is_dir():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    graph_params = yaml.safe_load(graph_params_path.read_text())
    calibrated_radius_um = graph_params["calibrated_radius_um"]
    graph_name = f"radius_{calibrated_radius_um}um"

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    final_annotations = pd.read_parquet(final_annotations_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    n_cross_patient_edges_total = 0

    for section_dir in sorted(d for d in pruned_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        graph_path = section_dir / f"{graph_name}.npz"
        if not graph_path.is_file():
            raise PipelineError(
                f"'{graph_path}' not found -- expected the `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`-calibrated graph to exist."
            )

        cell_order = pd.read_csv(section_dir / "cell_order.tsv", sep="\t")["cell_id"]
        coords = adata.obs.loc[cell_order, ["x_centroid", "y_centroid"]].to_numpy()
        adjacency = sparse.load_npz(graph_path)
        weighted = build_weighted_graph(adjacency, coords)

        node_metadata = pd.DataFrame(index=cell_order)
        node_metadata["patient_id"] = adata.obs.loc[cell_order, "patient_id"].to_numpy()
        node_metadata["section_id"] = section_id
        node_metadata["x_centroid"] = coords[:, 0]
        node_metadata["y_centroid"] = coords[:, 1]
        for col in NODE_METADATA_COLUMNS:
            node_metadata[col] = final_annotations.reindex(cell_order)[col].to_numpy()

        # Cross-check, not assumed: every edge must connect cells from
        # the same patient, since a section belongs to exactly one
        # patient by construction.
        coo = weighted.tocoo()
        patient_arr = node_metadata["patient_id"].to_numpy()
        n_cross_patient_edges = int((patient_arr[coo.row] != patient_arr[coo.col]).sum())
        n_cross_patient_edges_total += n_cross_patient_edges
        if n_cross_patient_edges > 0:
            raise PipelineError(
                f"'{section_id}': {n_cross_patient_edges} edge(s) connect different patients -- "
                "this violates the patient-separated-components invariant."
            )

        section_out_dir = output_dir / section_id
        section_out_dir.mkdir(parents=True, exist_ok=True)
        sparse.save_npz(section_out_dir / "primary_graph.npz", weighted)
        node_metadata.to_csv(section_out_dir / "node_metadata.tsv", sep="\t")

        manifest_rows.append(
            {
                "section_id": section_id,
                "patient_id": node_metadata["patient_id"].iloc[0],
                "n_nodes": weighted.shape[0],
                "n_edges": int(weighted.nnz / 2),
                "mean_degree": round(float(weighted.nnz / weighted.shape[0]), 4),
                "mean_edge_distance_um": round(float(coo.data.mean()), 4) if coo.nnz else None,
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    return {
        "graph_type": graph_name,
        "n_sections": len(manifest),
        "n_patients": manifest["patient_id"].nunique(),
        "n_nodes_total": int(manifest["n_nodes"].sum()),
        "n_edges_total": int(manifest["n_edges"].sum()),
        "n_cross_patient_edges_found": n_cross_patient_edges_total,
        "output_dir": str(output_dir),
        "manifest_path": str(manifest_path),
    }
