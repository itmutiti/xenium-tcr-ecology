"""Tumour-T-cell bipartite graph construction (`09_spatial_graph_construction_and_calibration/04_construct_tumour_tcell_bipartite_graph.py`).

Represents direct and near-direct malignant-cell/T-cell relationships for
contact-focused analyses (the concrete input `11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`, Clone Ecology Confirmatory Models's clone-
tumour engagement work needs), by extracting the cross-type edges already
present in `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`, `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`'s graphs rather than recomputing spatial
relationships from scratch.

**Malignant-cell population:** `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`'s `in_tumour_region` flag
(`data/derived/tumour_masks/*.parquet`) -- the most refined, spatially-
validated "this cell is part of a tumour region" call in this project
(built from `07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`'s malignancy score, spatially smoothed and
size-filtered), not the coarser `06_cell_type_annotation/06_integrate_annotation_evidence.py` `final_lineage ==
"Epithelial_Tumour"` alone. T-cell population: `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s
`final_lineage == "T_cell"`.

**Two contact tiers, matching the specification's "direct and near-
direct" phrasing:**
1. **Direct contact** -- the tumour-T-cell cross-type edges already
   present in `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`'s `boundary_contact` graph (touching cell-
   boundary polygons, 0.5um tolerance).
2. **Near-direct contact** -- the tumour-T-cell cross-type edges already
   present in `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated primary radius graph (30um,
   gap-pruned).

Both are extracted as subgraphs of the already-built, already-validated
graphs (not rebuilt from coordinates), keeping only edges with exactly
one malignant and one T-cell endpoint -- tumour-tumour, T-cell-T-cell,
and any edge touching a third cell type are all discarded, since this
graph exists specifically to represent the cross-type relationship.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError


def extract_bipartite_subgraph(
    graph: sparse.csr_matrix, is_malignant: np.ndarray, is_tcell: np.ndarray
) -> sparse.csr_matrix:
    """Pure, testable extraction of cross-type (malignant <-> T-cell)
    edges from an already-built graph. Same-type and third-type edges
    are dropped."""
    coo = graph.tocoo()
    cross = (is_malignant[coo.row] & is_tcell[coo.col]) | (
        is_tcell[coo.row] & is_malignant[coo.col]
    )
    return sparse.csr_matrix((coo.data[cross], (coo.row[cross], coo.col[cross])), shape=graph.shape)


def build_tumour_tcell_bipartite_graphs(project_root: Path) -> dict:
    pruned_dir = project_root / "data" / "graphs" / "pruned"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    tumour_masks_dir = project_root / "data" / "derived" / "tumour_masks"
    output_dir = project_root / "data" / "graphs" / "tumour_tcell"

    for p in (pruned_dir, final_annotations_path, tumour_masks_dir):
        if not p.exists():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    final_annotations = pd.read_parquet(final_annotations_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for section_dir in sorted(d for d in pruned_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        mask_path = tumour_masks_dir / f"{section_id}.parquet"
        if not mask_path.is_file():
            raise PipelineError(
                f"'{mask_path}' not found. Run `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py` first."
            )

        cell_order = pd.read_csv(section_dir / "cell_order.tsv", sep="\t")["cell_id"]
        tumour_mask_df = pd.read_parquet(mask_path)

        is_malignant = tumour_mask_df.reindex(cell_order)["in_tumour_region"].to_numpy(
            dtype=bool, na_value=False
        )
        lineage = final_annotations.reindex(cell_order)["final_lineage"]
        is_tcell = (lineage == "T_cell").to_numpy()

        section_out_dir = output_dir / section_id
        section_out_dir.mkdir(parents=True, exist_ok=True)
        (section_out_dir / "cell_order.tsv").write_text(
            (section_dir / "cell_order.tsv").read_text()
        )

        section_row = {
            "section_id": section_id,
            "n_malignant_cells": int(is_malignant.sum()),
            "n_tcells": int(is_tcell.sum()),
        }
        for tier, graph_name in [("direct", "boundary_contact"), ("near_direct", "radius_30.0um")]:
            graph = sparse.load_npz(section_dir / f"{graph_name}.npz")
            bipartite = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
            sparse.save_npz(section_out_dir / f"{tier}.npz", bipartite)
            section_row[f"n_{tier}_edges"] = int(bipartite.nnz / 2)
            section_row[f"n_tcells_with_{tier}_contact"] = int(
                ((np.asarray(bipartite.sum(axis=1)).ravel() > 0) & is_tcell).sum()
            )

        manifest_rows.append(section_row)

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_sections": len(manifest),
        "n_malignant_cells_total": int(manifest["n_malignant_cells"].sum()),
        "n_tcells_total": int(manifest["n_tcells"].sum()),
        "n_direct_contact_edges_total": int(manifest["n_direct_edges"].sum()),
        "n_near_direct_contact_edges_total": int(manifest["n_near_direct_edges"].sum()),
        "n_tcells_with_any_direct_contact": int(manifest["n_tcells_with_direct_contact"].sum()),
        "n_tcells_with_any_near_direct_contact": int(
            manifest["n_tcells_with_near_direct_contact"].sum()
        ),
        "output_dir": str(output_dir),
        "manifest_path": str(manifest_path),
    }
