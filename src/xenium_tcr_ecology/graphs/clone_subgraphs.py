"""Clone-induced subgraph construction (`09_spatial_graph_construction_and_calibration/05_construct_clone_induced_subgraphs.py`).

Extracts all cells belonging to each high-confidence clone (`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`'s
frozen `high_confidence_clones.parquet`) plus local microenvironment
shells -- concentric graph-distance rings around the clone's own cells in
`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated primary graph -- the concrete spatial-context
input Niche and Ecosystem Discovery+ niche/microenvironment analyses need per clone.

**Clone membership:** re-derived from `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`'s
`tcr_resolved_calls.parquet` (`detected_probes` == the frozen clone's
`clone_id`, restricted to `singlet`/`low_confidence` resolution, matching
`08_tcr_clonal_analysis/07_build_clone_metadata_table.py`'s own clone definition exactly) -- `clone_metadata.parquet`
itself only stores clone-level aggregates, not the member cell list, so
membership must be reconstructed from its original source, not
approximated.

**Microenvironment shells:** `N_SHELLS = 3` concentric graph-distance
rings (BFS expansion over `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated 30um primary graph, not
a fixed physical radius -- shell 1 is "directly graph-connected to a
clone member," shell 2 is "two graph hops away," etc., which naturally
adapts to local tissue density the same way the calibrated primary graph
itself does). Shell 0 is the clone's own member cells. A cell already in
the clone (shell 0) or a closer shell is never re-labelled into a farther
shell if reached by a longer path from a different clone member (standard
BFS shortest-hop-distance semantics).

Clones spanning multiple sections (technical-replicate patients,
`08_tcr_clonal_analysis/07_build_clone_metadata_table.py`) are expanded independently per section -- `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s
primary graph never connects cells across sections, so a clone's shells
in one section cannot depend on its cells in another.
"""

from __future__ import annotations

from collections import deque
from pathlib import Path

import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

N_SHELLS = 3


def compute_clone_shells(
    clone_cell_indices: set[int], graph: sparse.csr_matrix, n_shells: int = N_SHELLS
) -> dict[int, int]:
    """Pure, testable BFS shell expansion. Returns {cell_index:
    shell_distance} for the clone's own cells (0) and every cell reached
    within `n_shells` graph hops (1..n_shells) -- cells farther than
    `n_shells` hops, or in a disconnected component, are absent."""
    distances: dict[int, int] = {idx: 0 for idx in clone_cell_indices}
    frontier = deque(clone_cell_indices)
    while frontier:
        node = frontier.popleft()
        current_distance = distances[node]
        if current_distance >= n_shells:
            continue
        neighbors = graph.indices[graph.indptr[node] : graph.indptr[node + 1]]
        for neighbor in neighbors:
            if neighbor not in distances:
                distances[neighbor] = current_distance + 1
                frontier.append(neighbor)
    return distances


def build_clone_subgraphs(project_root: Path) -> dict:
    release_dir = project_root / "data" / "releases" / "v1_tcr_calls"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    output_dir = project_root / "data" / "graphs" / "clones"

    high_confidence_path = release_dir / "high_confidence_clones.parquet"
    resolved_calls_path = release_dir / "tcr_resolved_calls.parquet"
    for p in (high_confidence_path, resolved_calls_path, primary_graphs_dir):
        if not p.exists():
            raise PipelineError(
                f"'{p}' not found. Run `08_tcr_clonal_analysis/08_generate_tcr_release_report.py` and `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` first."
            )

    high_confidence = pd.read_parquet(high_confidence_path)
    resolved_calls = pd.read_parquet(resolved_calls_path)
    clonal_cells = resolved_calls[resolved_calls["resolution"].isin(["singlet", "low_confidence"])]

    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, clone_row in high_confidence.iterrows():
        clone_id = clone_row["clone_id"]
        member_cell_ids = clonal_cells.index[clonal_cells["detected_probes"] == clone_id]
        if len(member_cell_ids) == 0:
            raise PipelineError(
                f"Clone '{clone_id}' has 0 member cells in tcr_resolved_calls.parquet -- data inconsistency."
            )

        by_section = resolved_calls.loc[member_cell_ids].groupby("section_id", observed=True)
        for section_id, section_group in by_section:
            section_dir = primary_graphs_dir / section_id
            node_metadata = pd.read_csv(section_dir / "node_metadata.tsv", sep="\t", index_col=0)
            graph = sparse.load_npz(section_dir / "primary_graph.npz")

            cell_id_to_index = {cid: i for i, cid in enumerate(node_metadata.index)}
            clone_indices = {
                cell_id_to_index[cid] for cid in section_group.index if cid in cell_id_to_index
            }
            if not clone_indices:
                continue

            shells = compute_clone_shells(clone_indices, graph)
            index_to_cell_id = {i: cid for cid, i in cell_id_to_index.items()}
            shell_df = pd.DataFrame(
                {
                    "cell_id": [index_to_cell_id[i] for i in shells],
                    "shell": list(shells.values()),
                }
            ).set_index("cell_id")

            safe_clone_id = clone_id.replace(";", "__").replace("/", "_")
            clone_section_dir = output_dir / safe_clone_id
            clone_section_dir.mkdir(parents=True, exist_ok=True)
            shell_df.to_csv(clone_section_dir / f"{section_id}.tsv", sep="\t")

            rows.append(
                {
                    "clone_id": clone_id,
                    "section_id": section_id,
                    "n_clone_cells": int((shell_df["shell"] == 0).sum()),
                    **{
                        f"n_shell{s}_cells": int((shell_df["shell"] == s).sum())
                        for s in range(1, N_SHELLS + 1)
                    },
                }
            )

    manifest = pd.DataFrame(rows)
    manifest_path = output_dir / "_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_clones": manifest["clone_id"].nunique(),
        "n_clone_section_instances": len(manifest),
        "median_n_clone_cells": float(manifest["n_clone_cells"].median()),
        "median_n_shell1_cells": float(manifest["n_shell1_cells"].median()),
        "output_dir": str(output_dir),
        "manifest_path": str(manifest_path),
    }
