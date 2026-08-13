"""Graph parameter calibration (`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`).

Selects one calibrated radius and one calibrated k (from `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`'s
candidates) for `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s "primary" cell graph, using
connectivity evidence measured on the `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py` pruned graphs across all
18 sections -- not asserted qualitatively.

**Evidence, computed before choosing anything:**
connected-component analysis (`scipy.sparse.csgraph.connected_components`,
the same tool already used for tumour-region labelling in `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`) on
every candidate graph, every section. Result: `radius_15um` is badly
fragmented (median largest-component fraction 0.836, as low as 0.121 in
one section) despite being the most gap-robust candidate (Phase 9.01:
0% pruned) -- it is simply too sparse on its own to form a usable
connected graph. `radius_30um` (median 0.985) and `radius_50um` (median
0.988) are both well-connected post-pruning; all three k-NN candidates
are essentially always fully connected by construction (median >=0.999).

**Selection rule:** `CONNECTIVITY_THRESHOLD = 0.95` (a standard,
commonly-used bar for "a graph is effectively connected" -- most cells in
one giant component, not scattered across many small ones). Among radius
candidates clearing this bar, the smallest is preferred (more locally
interpretable as a genuine short-range interaction scale, and Phase
9.01's evidence showed larger radii are progressively more contaminated
by gap-bridging edges before pruning, even though pruning itself resolves
this): `radius_30um` is selected over `radius_50um`. Among k-NN
candidates (all fully connected), the smallest, k=6, is selected for the
same "most immediate, most locally interpretable" reasoning --
consistent with `RADIUS_CANDIDATES_UM`'s data grounding (`09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`) in the
median cell diameter (9.82um): a 6-neighbour k-NN graph approximates one
ring of immediate physical neighbours, the same biological scale
radius_15um was intended to capture but could not deliver on its own due
to fragmentation.

**Finding that does not change the calibration decision:** two sections
(P10_run1, P19_run2) show genuine multi-fragment structure even in the
fully-connected k-NN graphs (P10_run1: two comparable-sized components of
19,571 and 18,351 cells, roughly 52%/48%, plus one 7-cell fragment) --
physically separate tissue pieces within one section, not a
graph-construction defect (confirmed because even k-NN, which is
essentially guaranteed connected for a single contiguous point cloud,
shows the same split). This is expected, documented pipeline behaviour,
not a blocker: `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s spec already anticipates
"patient-separated components."
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from xenium_tcr_ecology.infra.exceptions import PipelineError

CONNECTIVITY_THRESHOLD = 0.95


def compute_largest_component_fraction(graph: sparse.csr_matrix) -> float:
    """Pure, testable connectivity summary for one graph."""
    if graph.shape[0] == 0:
        return 1.0
    _, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    return float(sizes.max() / graph.shape[0])


def select_calibrated_parameter(
    connectivity_by_candidate: dict[str, float],
    candidate_scale: dict[str, float],
    threshold: float = CONNECTIVITY_THRESHOLD,
) -> str:
    """Pure, testable selection rule: among candidates whose median
    connectivity clears `threshold`, pick the one with the smallest
    underlying scale (`candidate_scale`, e.g. radius in um or k); if none
    clear the bar, fall back to the single most-connected candidate."""
    passing = [c for c, frac in connectivity_by_candidate.items() if frac >= threshold]
    if passing:
        return min(passing, key=lambda c: candidate_scale[c])
    return max(connectivity_by_candidate, key=connectivity_by_candidate.get)


def build_graph_parameter_calibration(project_root: Path) -> dict:
    pruned_dir = project_root / "data" / "graphs" / "pruned"
    output_path = project_root / "config" / "graph_parameters.yaml"

    if not pruned_dir.is_dir():
        raise PipelineError(
            f"'{pruned_dir}' not found. Run `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py` first."
        )

    radius_candidates = {"radius_15.0um": 15.0, "radius_30.0um": 30.0, "radius_50.0um": 50.0}
    knn_candidates = {"knn_6": 6, "knn_10": 10, "knn_15": 15}
    all_candidates = {**radius_candidates, **knn_candidates}

    section_dirs = sorted(d for d in pruned_dir.iterdir() if d.is_dir())
    if not section_dirs:
        raise PipelineError(f"No section directories found under '{pruned_dir}'.")

    rows = []
    for section_dir in section_dirs:
        for graph_name in all_candidates:
            graph_path = section_dir / f"{graph_name}.npz"
            if not graph_path.is_file():
                continue
            graph = sparse.load_npz(graph_path)
            frac = compute_largest_component_fraction(graph)
            rows.append(
                {
                    "section_id": section_dir.name,
                    "graph_name": graph_name,
                    "largest_component_fraction": frac,
                }
            )

    connectivity_df = pd.DataFrame(rows)
    median_connectivity = (
        connectivity_df.groupby("graph_name")["largest_component_fraction"].median().to_dict()
    )

    radius_connectivity = {k: v for k, v in median_connectivity.items() if k in radius_candidates}
    knn_connectivity = {k: v for k, v in median_connectivity.items() if k in knn_candidates}

    calibrated_radius_name = select_calibrated_parameter(radius_connectivity, radius_candidates)
    calibrated_knn_name = select_calibrated_parameter(knn_connectivity, knn_candidates)

    # Genuine multi-fragment tissue structure is only meaningfully
    # signalled by a k-NN graph fragmenting: k-NN connects every cell to
    # its k nearest neighbours regardless of absolute distance, so it is
    # essentially guaranteed connected for any single contiguous point
    # cloud -- confirmed: an earlier version of this
    # check spanned all six candidate graph types, which mostly just
    # re-flagged radius_15um's already-known, uninteresting sparsity
    # (14/18 sections), not tissue separation. Restricting to k-NN
    # candidates only gives the intended, specific signal.
    knn_connectivity_df = connectivity_df[connectivity_df["graph_name"].isin(knn_candidates)]
    fragmented_sections = (
        knn_connectivity_df[knn_connectivity_df["largest_component_fraction"] < 0.9]["section_id"]
        .unique()
        .tolist()
    )

    config = {
        "calibrated_radius_um": radius_candidates[calibrated_radius_name],
        "calibrated_knn": knn_candidates[calibrated_knn_name],
        "connectivity_threshold": CONNECTIVITY_THRESHOLD,
        "median_largest_component_fraction_by_candidate": {
            k: round(v, 4) for k, v in median_connectivity.items()
        },
        "sections_with_genuine_multi_fragment_structure": sorted(fragmented_sections),
        "rationale": (
            "Selected the smallest radius/k whose median largest-component fraction across all 18 sections "
            "clears 0.95, preferring locally interpretable short-range scales over larger ones once both are "
            "adequately connected."
        ),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(config, sort_keys=False))

    return {
        "calibrated_radius_um": config["calibrated_radius_um"],
        "calibrated_knn": config["calibrated_knn"],
        "median_largest_component_fraction_by_candidate": config[
            "median_largest_component_fraction_by_candidate"
        ],
        "n_sections_with_genuine_multi_fragment_structure": len(fragmented_sections),
        "output_path": str(output_path),
    }
