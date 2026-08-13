"""Clone-specific stromal/myeloid barrier metrics (`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`).

For each (clone, section) -- `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`'s established unit -- barrier
composition along the shortest graph path from a clone's cells to the
nearest `Epithelial_Tumour` cell, for three barrier lineage groups:

- **fibroblast**: the full `Fibroblast` lineage (`Activated_CAF` +
  `Resting_Fibroblast` + unresolved substate).
- **vascular**: `Endothelial` + `Perivascular_SmoothMuscle` lineages --
  the same pairing `10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py` found strongly co-enriched (z=39.96).
- **suppressive_myeloid**: the `Macrophage` myeloid substate
  specifically. This project's `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R` taxonomy has no further
  M1/M2 polarisation resolution -- `Macrophage` is the literature-
  standard "suppressive myeloid" proxy this project's substate taxonomy
  supports; `Monocyte` is deliberately excluded (a weaker, less
  literature-established "suppressive" proxy) rather than fabricating a
  finer polarisation category.

**Method:** a single multi-source Dijkstra shortest-path solve per
section (`scipy.sparse.csgraph.dijkstra`, `min_only=True`,
`return_predecessors=True`) from every `Epithelial_Tumour` cell
simultaneously, on `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated, distance-weighted (micron
edge weights, confirmed) primary graph -- not one solve per
clone cell. Every cell's predecessor chain traces its shortest path to
whichever tumour cell is nearest; a clone cell's intermediate path cells
(excluding both endpoints) give the barrier composition for that cell.

**Calibrated null model, reusing `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s established null-model
type:** constrained permutation of lineage labels within a section (the
same conceptual null `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py` already calibrated for Type I error/power
on synthetic ground truth) -- not re-solved per permutation. Because
shortest paths depend only on graph topology and distance weights, never
on lineage labels, the paths are traced once per section and reused
across all `N_PERMUTATIONS = 199` draws (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s established
count); only the lineage-label assignment is reshuffled per permutation
(a cheap array operation), giving an efficient empirical null without
re-solving Dijkstra 199 times.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra

from xenium_tcr_ecology.clone_ecology.spatial_descriptors import filter_to_primary_cohort
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed

FIBROBLAST_LINEAGES = {"Fibroblast"}
VASCULAR_LINEAGES = {"Endothelial", "Perivascular_SmoothMuscle"}
SUPPRESSIVE_MYELOID_SUBSTATE = "Macrophage"
TUMOUR_LINEAGE = "Epithelial_Tumour"
BARRIER_GROUPS = ["fibroblast", "vascular", "suppressive_myeloid"]
N_PERMUTATIONS = 199
RNG_SEED = get_default_seed()
UNREACHABLE_PREDECESSOR = -9999


def assign_barrier_group(lineage: np.ndarray, substate: np.ndarray) -> np.ndarray:
    """Pure, testable: per-cell barrier-group label ('fibroblast',
    'vascular', 'suppressive_myeloid', or 'other') from lineage and
    substate arrays."""
    group = np.full(len(lineage), "other", dtype=object)
    group[np.isin(lineage, list(FIBROBLAST_LINEAGES))] = "fibroblast"
    group[np.isin(lineage, list(VASCULAR_LINEAGES))] = "vascular"
    group[substate == SUPPRESSIVE_MYELOID_SUBSTATE] = "suppressive_myeloid"
    return group


def trace_intermediate_path(
    predecessors: np.ndarray, target_idx: int, source_idx_set: set[int]
) -> list[int]:
    """Pure, testable: intermediate node indices on the shortest path
    from `target_idx` back to whichever source (tumour cell) achieved
    its shortest distance, excluding both `target_idx` itself and the
    terminal source cell. Empty list if unreachable."""
    path: list[int] = []
    current = target_idx
    while True:
        prev = predecessors[current]
        if prev == UNREACHABLE_PREDECESSOR:
            return []
        if prev in source_idx_set:
            return path
        path.append(prev)
        current = prev


def compute_barrier_fraction(path_barrier_groups: np.ndarray, barrier_group: str) -> float:
    """Pure, testable: fraction of intermediate path cells belonging to
    `barrier_group`. NaN for an empty (unreachable) path."""
    if len(path_barrier_groups) == 0:
        return float("nan")
    return float(np.mean(path_barrier_groups == barrier_group))


def compute_empirical_pvalue_greater(observed: float, null_values: np.ndarray) -> float:
    """Pure, testable: one-sided empirical p-value for observed barrier
    fraction being unusually high relative to the null (label-permuted)
    distribution -- same convention as `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s calibrated
    permutation tests."""
    n = len(null_values)
    return float((np.sum(np.asarray(null_values) >= observed) + 1) / (n + 1))


def build_clone_barrier_metrics(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    output_path = project_root / "data" / "derived" / "clone_barrier_metrics.parquet"

    for path, phase in [
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
        (sample_manifest_path, None),
        (final_annotations_path, "`06_cell_type_annotation/06_integrate_annotation_evidence.py`"),
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
            "No primary-cohort, release clone cells found -- cannot build barrier metrics."
        )

    final_annotations = pd.read_parquet(final_annotations_path)
    rng = np.random.default_rng(RNG_SEED)

    rows = []
    for section_id in sorted(clonal["section_id"].unique()):
        node_metadata = pd.read_csv(primary_graphs_dir / section_id / "node_metadata.tsv", sep="\t")
        graph = sparse.load_npz(primary_graphs_dir / section_id / "primary_graph.npz")

        section_annotations = final_annotations.reindex(node_metadata["cell_id"])
        lineage = section_annotations["final_lineage"].to_numpy()
        substate = section_annotations["final_substate"].to_numpy()
        real_barrier_group = assign_barrier_group(lineage, substate)

        tumour_idx = np.flatnonzero(lineage == TUMOUR_LINEAGE)
        if len(tumour_idx) == 0:
            continue
        tumour_idx_set = set(tumour_idx.tolist())

        # min_only=True + return_predecessors=True returns a 3-tuple
        # (dist, predecessors, sources) -- confirmed against the
        # installed scipy version before relying on it (a naive 2-tuple
        # unpack raises ValueError here).
        _, predecessors, _ = dijkstra(
            graph, directed=False, indices=tumour_idx, min_only=True, return_predecessors=True
        )

        cell_id_to_idx = {cid: i for i, cid in enumerate(node_metadata["cell_id"])}
        section_clonal = clonal[clonal["section_id"] == section_id]

        # Paths are computed once per section (predecessors above); the
        # null permutes barrier-group labels only (cheap) and reuses the
        # same path indices for every permutation.
        null_barrier_groups = [rng.permutation(real_barrier_group) for _ in range(N_PERMUTATIONS)]

        for clone_id, group in section_clonal.groupby("clone_id", observed=True):
            cell_indices = [cell_id_to_idx[cid] for cid in group.index if cid in cell_id_to_idx]
            if len(cell_indices) == 0:
                continue

            observed_fractions: dict[str, list[float]] = {
                barrier_group: [] for barrier_group in BARRIER_GROUPS
            }
            null_fractions = {
                barrier_group: np.zeros(N_PERMUTATIONS) for barrier_group in BARRIER_GROUPS
            }

            all_intermediate_idx: list[int] = []
            for cell_idx in cell_indices:
                path_idx = trace_intermediate_path(predecessors, cell_idx, tumour_idx_set)
                all_intermediate_idx.extend(path_idx)
                path_groups = real_barrier_group[path_idx] if path_idx else np.array([])
                for barrier_group in BARRIER_GROUPS:
                    observed_fractions[barrier_group].append(
                        compute_barrier_fraction(path_groups, barrier_group)
                    )

            all_intermediate_idx_arr = np.array(all_intermediate_idx, dtype=int)
            # A clone whose cells are all unreachable to any tumour cell
            # (a disconnected-component possibility, `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`) has an
            # all-NaN observed_fractions list -- checked explicitly
            # rather than relying on np.nanmean's silent-but-warning
            # all-NaN reduction (same fix as `10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py`).
            mean_observed = {}
            for barrier_group in BARRIER_GROUPS:
                values = np.asarray(observed_fractions[barrier_group], dtype=float)
                mean_observed[barrier_group] = (
                    float(np.nanmean(values)) if np.isfinite(values).any() else float("nan")
                )

            if len(all_intermediate_idx_arr) > 0:
                for p in range(N_PERMUTATIONS):
                    perm_groups = null_barrier_groups[p][all_intermediate_idx_arr]
                    for barrier_group in BARRIER_GROUPS:
                        null_fractions[barrier_group][p] = np.mean(perm_groups == barrier_group)
                pvalues = {
                    barrier_group: compute_empirical_pvalue_greater(
                        mean_observed[barrier_group], null_fractions[barrier_group]
                    )
                    for barrier_group in BARRIER_GROUPS
                }
            else:
                pvalues = {barrier_group: float("nan") for barrier_group in BARRIER_GROUPS}

            rows.append(
                {
                    "clone_id": clone_id,
                    "section_id": section_id,
                    "patient_id": group["patient_id"].iloc[0],
                    "n_cells": len(group),
                    "n_cells_with_path": len(cell_indices),
                    "fibroblast_barrier_fraction": mean_observed["fibroblast"],
                    "fibroblast_barrier_pvalue": pvalues["fibroblast"],
                    "vascular_barrier_fraction": mean_observed["vascular"],
                    "vascular_barrier_pvalue": pvalues["vascular"],
                    "suppressive_myeloid_barrier_fraction": mean_observed["suppressive_myeloid"],
                    "suppressive_myeloid_barrier_pvalue": pvalues["suppressive_myeloid"],
                }
            )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_clone_section_rows": len(result),
        "n_distinct_clones": int(result["clone_id"].nunique()) if len(result) else 0,
        "n_significant_fibroblast_barrier": (
            int((result["fibroblast_barrier_pvalue"] < 0.05).sum()) if len(result) else 0
        ),
        "n_significant_vascular_barrier": (
            int((result["vascular_barrier_pvalue"] < 0.05).sum()) if len(result) else 0
        ),
        "n_significant_suppressive_myeloid_barrier": (
            int((result["suppressive_myeloid_barrier_pvalue"] < 0.05).sum()) if len(result) else 0
        ),
        "output_path": str(output_path),
    }
