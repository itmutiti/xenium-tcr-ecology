"""Tissue domain segmentation (`10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`).

Forms spatially contiguous tissue domains from `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s per-cell
archetype labels, while preserving genuine boundaries and rare niches --
not a from-scratch spatial clustering, but a principled cleanup of the
raw per-cell k-means archetype call.

**Method:** majority-vote smoothing of archetype labels over `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s
calibrated 30um primary graph -- the same scale already used for
archetype discovery itself (`10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`), not a new, independently-chosen
smoothing neighbourhood. Because each cell's archetype label already
reflects its 30um neighbourhood composition (`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`), the raw label
field is already spatially coherent; majority-vote smoothing over that
same local neighbourhood only removes residual per-cell instability from
the discrete k-means assignment (isolated "salt-and-pepper" mislabels),
not genuine spatial structure -- so boundaries are not blurred. Ties are
broken in favour of the cell's own original label (a defensible,
deterministic rule, avoiding arbitrary numeric-label-ordering bias at
genuine near-50/50 boundary cells).

A tissue domain is then a maximal connected set of same-smoothed-label
cells under that same graph -- connected components of the graph
restricted to same-label edges. Two spatially separate patches sharing a
label become two distinct domains, never merged. No minimum-size
threshold is applied anywhere: an isolated single-cell or
few-cell patch keeps its own `domain_id`, preserving rare niches rather
than merging or discarding them.

The 61 (all-sections total) cells with zero degree in `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s
primary graph have no `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R` archetype label (`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`'s
zero-degree exclusion at this same scale) and are excluded here too, for
the same reason: "no neighbours to describe" is a different
condition from "neighbours exist," not one to paper over with a
fabricated label.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from xenium_tcr_ecology.infra.exceptions import PipelineError


def smooth_labels_by_majority_vote(graph: sparse.csr_matrix, labels: np.ndarray) -> np.ndarray:
    """Pure, testable per-cell majority vote over {self} union graph
    neighbours. `labels` must be integer-coded (0..k-1). Ties (including
    the case where the cell's own label already achieves the max vote
    count) are broken in favour of the cell's own original label."""
    n = graph.shape[0]
    k = int(labels.max()) + 1

    graph_binary = graph.copy()
    graph_binary.data = np.ones_like(graph_binary.data)

    one_hot = np.zeros((n, k), dtype=np.float64)
    one_hot[np.arange(n), labels] = 1.0
    total_votes = (graph_binary @ one_hot) + one_hot  # neighbours + self

    own_label_votes = total_votes[np.arange(n), labels]
    max_votes = total_votes.max(axis=1)
    own_label_is_a_mode = own_label_votes >= max_votes - 1e-9
    majority_label = total_votes.argmax(axis=1)

    return np.where(own_label_is_a_mode, labels, majority_label)


def find_contiguous_domains(graph: sparse.csr_matrix, labels: np.ndarray) -> np.ndarray:
    """Pure, testable: connected components of `graph` restricted to
    edges connecting cells with the same label. Domains are never merged
    across disconnected same-label patches; no size floor is applied."""
    coo = graph.tocoo()
    same_label = labels[coo.row] == labels[coo.col]
    filtered = sparse.coo_matrix(
        (coo.data[same_label], (coo.row[same_label], coo.col[same_label])),
        shape=graph.shape,
    ).tocsr()
    _, domain_ids = connected_components(filtered, directed=False)
    return domain_ids


def build_tissue_domains(project_root: Path) -> dict:
    archetypes_path = project_root / "data" / "derived" / "neighbourhood_archetypes.parquet"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    output_path = project_root / "data" / "derived" / "tissue_domains.parquet"

    if not archetypes_path.is_file():
        raise PipelineError(
            f"'{archetypes_path}' not found. Run `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R` first."
        )
    if not primary_graphs_dir.is_dir():
        raise PipelineError(
            f"'{primary_graphs_dir}' not found. Run `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` first."
        )

    archetypes = pd.read_parquet(archetypes_path).set_index("cell_id")
    archetype_labels_order = sorted(archetypes["archetype"].unique())
    label_to_code = {label: i for i, label in enumerate(archetype_labels_order)}
    code_to_label = {i: label for label, i in label_to_code.items()}

    section_results = []
    n_cells_excluded_zero_degree = 0
    n_cells_relabelled_by_smoothing = 0
    domain_id_offset = 0
    for section_dir in sorted(d for d in primary_graphs_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        node_metadata_path = section_dir / "node_metadata.tsv"
        graph_path = section_dir / "primary_graph.npz"
        if not node_metadata_path.is_file() or not graph_path.is_file():
            raise PipelineError(
                f"'{section_dir}' is missing 'node_metadata.tsv' or 'primary_graph.npz'. Run `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` first."
            )

        node_metadata = pd.read_csv(node_metadata_path, sep="\t")
        graph = sparse.load_npz(graph_path)

        section_archetype = archetypes.reindex(node_metadata["cell_id"])["archetype"]
        included_mask = section_archetype.notna().to_numpy()
        n_cells_excluded_zero_degree += int((~included_mask).sum())

        included_idx = np.flatnonzero(included_mask)
        sub_graph = graph[included_idx][:, included_idx]
        raw_codes = section_archetype.to_numpy()[included_idx].astype(int)
        raw_codes = np.array([label_to_code[v] for v in raw_codes])

        smoothed_codes = smooth_labels_by_majority_vote(sub_graph, raw_codes)
        n_cells_relabelled_by_smoothing += int((smoothed_codes != raw_codes).sum())

        local_domain_ids = find_contiguous_domains(sub_graph, smoothed_codes)
        global_domain_ids = local_domain_ids + domain_id_offset
        domain_id_offset += int(local_domain_ids.max()) + 1 if len(local_domain_ids) > 0 else 0

        section_result = pd.DataFrame(
            {
                "cell_id": node_metadata["cell_id"].to_numpy()[included_idx],
                "section_id": section_id,
                "archetype": [code_to_label[c] for c in smoothed_codes],
                "domain_id": global_domain_ids,
            }
        )
        section_results.append(section_result)

    result = pd.concat(section_results, ignore_index=True)
    domain_sizes = result.groupby("domain_id")["cell_id"].transform("size")
    result["domain_size"] = domain_sizes

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_cells": len(result),
        "n_cells_excluded_zero_degree": n_cells_excluded_zero_degree,
        "n_cells_relabelled_by_smoothing": n_cells_relabelled_by_smoothing,
        "n_domains": int(result["domain_id"].nunique()),
        "n_single_cell_domains": int((domain_sizes == 1).sum()),
        "output_path": str(output_path),
    }
