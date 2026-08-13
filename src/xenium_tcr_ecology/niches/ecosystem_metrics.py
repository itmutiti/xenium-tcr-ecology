"""Ecosystem abundance and topology metrics (`10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py`).

Quantifies six distinct spatial properties of each `10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`
ecosystem, per section, per the scaffold's specification (abundance,
fragmentation, interface length, mixing, compactness, relation to
tumour borders) -- each grounded directly in an already-computed
artifact from earlier phases, not a new, independently invented geometry
pipeline:

1. **Abundance**: fraction of a section's cells belonging to the
   ecosystem.
2. **Fragmentation**: `n_domains` (`10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`'s domain count for this
   ecosystem in this section) and `effective_n_patches` -- the inverse
   Simpson index of domain sizes (a standard landscape-ecology
   "effective number of equally-sized patches" statistic), which
   corrects raw domain count for the case where one dominant domain plus
   many tiny fragments would otherwise look identical to many
   similarly-sized domains.
3. **Interface length**: count of `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` primary-graph edges
   connecting a cell of this ecosystem to a cell of a different
   ecosystem within the section -- a standard, defensible proxy for
   interface length in a roughly uniform-density spatial point pattern
   (edge count scales with shared boundary length), not a literal
   polygon perimeter computation this dataset has no cell-shape polygons
   to support.
4. **Mixing**: mean fraction of this ecosystem's cells' graph
   neighbours that belong to a different ecosystem, averaged over all of
   this ecosystem's cells in the section (0 = fully segregated/interior,
   1 = every neighbour is a different ecosystem). Uses the same
   binarized `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` graph as the interface-length metric, not a new
   neighbourhood definition.
5. **Compactness**: cell-count-weighted mean isoperimetric quotient
   (4*pi*Area/Perimeter^2 of each domain's convex hull, from cell
   centroid coordinates in `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s `node_metadata.tsv`) across this
   ecosystem's domains with >=3 cells in the section (a value of 1 is a
   perfect circle; lower values are more elongated/irregular). Domains
   with <3 cells cannot form a hull and are excluded from the weighted
   mean, not assigned a fabricated value.
6. **Relation to tumour borders**: mean and median `signed_distance_to_
   tumour_boundary_um`, plus `fraction_inside_tumour_region`, all joined
   against Tumour Epithelium Characterisation's `data/derived/tumour_boundaries.parquet` by cell
   ID -- reused directly, not recomputed. Both mean and median are
   reported deliberately: this field is heavily right-skewed (Tumour Epithelium Characterisation's
   tumour mask is a tight, point-buffered union that "hugs" individual
   cell clusters rather than approximating one smoothed macro-scale
   margin; a small fraction of far-
   outlying cells can inflate the mean well above the typical cell's
   distance), so the mean alone would be a materially misleading summary.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import ConvexHull, QhullError

from xenium_tcr_ecology.infra.exceptions import PipelineError


def compute_effective_n_patches(domain_sizes: np.ndarray) -> float:
    """Pure, testable: inverse Simpson index of domain sizes."""
    total = domain_sizes.sum()
    if total == 0:
        return 0.0
    proportions = domain_sizes / total
    return float(1.0 / np.sum(proportions**2))


def compute_interface_edge_count(
    graph_binary: sparse.csr_matrix, ecosystem_codes: np.ndarray, this_code: int
) -> int:
    """Pure, testable: count of undirected graph edges (each counted
    once) connecting a `this_code` cell to a different-code cell."""
    coo = graph_binary.tocoo()
    upper = coo.row < coo.col
    row, col = coo.row[upper], coo.col[upper]
    is_this_row = ecosystem_codes[row] == this_code
    is_this_col = ecosystem_codes[col] == this_code
    crosses = is_this_row != is_this_col  # exactly one endpoint is `this_code`
    return int(np.sum(crosses))


def compute_mixing_index(
    graph_binary: sparse.csr_matrix, ecosystem_codes: np.ndarray, this_code: int
) -> float:
    """Pure, testable: mean fraction of a this-ecosystem cell's graph
    neighbours that are a different ecosystem, averaged over all
    this-ecosystem cells with nonzero degree."""
    mask = ecosystem_codes == this_code
    if not mask.any():
        return float("nan")

    sub = graph_binary[mask]
    degree = np.asarray(sub.sum(axis=1)).ravel()
    same_code_indicator = (ecosystem_codes == this_code).astype(np.float64)
    same_neighbour_count = np.asarray(sub @ same_code_indicator).ravel()
    diff_neighbour_count = degree - same_neighbour_count

    valid = degree > 0
    if not valid.any():
        return float("nan")
    with np.errstate(invalid="ignore", divide="ignore"):
        mixing_per_cell = diff_neighbour_count[valid] / degree[valid]
    return float(np.mean(mixing_per_cell))


def compute_isoperimetric_quotient(x: np.ndarray, y: np.ndarray) -> float | None:
    """Pure, testable: 4*pi*Area/Perimeter^2 of the convex hull of
    (x, y) points. Returns None for <3 points or a degenerate
    (collinear) hull -- never a fabricated value."""
    if len(x) < 3:
        return None
    points = np.column_stack([x, y])
    try:
        hull = ConvexHull(points)
    except QhullError:
        return None
    area = hull.volume  # scipy 2D convention: `volume` is area
    perimeter = hull.area  # scipy 2D convention: `area` is perimeter
    if perimeter <= 0:
        return None
    return float(4 * np.pi * area / perimeter**2)


def build_ecosystem_metrics(project_root: Path) -> dict:
    domains_path = project_root / "data" / "derived" / "tissue_domains.parquet"
    annotation_path = project_root / "metadata" / "ecosystem_annotation.tsv"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    tumour_boundaries_path = project_root / "data" / "derived" / "tumour_boundaries.parquet"
    output_path = project_root / "data" / "derived" / "ecosystem_metrics.parquet"

    for path, phase in [
        (domains_path, "`10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`"),
        (
            annotation_path,
            "`10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`",
        ),
        (
            primary_graphs_dir,
            "`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`",
        ),
        (tumour_boundaries_path, "Tumour Epithelium Characterisation"),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    domains = pd.read_parquet(domains_path).set_index("cell_id")
    annotation = pd.read_csv(annotation_path, sep="\t").set_index("archetype")["ecosystem_label"]
    domains["ecosystem_label"] = domains["archetype"].map(annotation)
    tumour_boundaries = pd.read_parquet(tumour_boundaries_path)

    ecosystem_labels_order = sorted(annotation.unique())
    label_to_code = {label: i for i, label in enumerate(ecosystem_labels_order)}

    section_ecosystem_rows = []
    for section_dir in sorted(d for d in primary_graphs_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        node_metadata = pd.read_csv(section_dir / "node_metadata.tsv", sep="\t")
        graph = sparse.load_npz(section_dir / "primary_graph.npz")
        graph_binary = graph.copy()
        graph_binary.data = np.ones_like(graph_binary.data)

        section_domains = domains.reindex(node_metadata["cell_id"])
        included_mask = section_domains["ecosystem_label"].notna().to_numpy()
        included_idx = np.flatnonzero(included_mask)
        sub_graph_binary = graph_binary[included_idx][:, included_idx]

        cell_id = node_metadata["cell_id"].to_numpy()[included_idx]
        x = node_metadata["x_centroid"].to_numpy()[included_idx]
        y = node_metadata["y_centroid"].to_numpy()[included_idx]
        ecosystem_label = section_domains["ecosystem_label"].to_numpy()[included_idx]
        ecosystem_codes = np.array([label_to_code[v] for v in ecosystem_label])
        domain_id = section_domains["domain_id"].to_numpy()[included_idx]
        domain_size = section_domains["domain_size"].to_numpy()[included_idx]

        section_tumour = tumour_boundaries.reindex(cell_id)
        n_cells_section = len(cell_id)

        for label in ecosystem_labels_order:
            code = label_to_code[label]
            mask = ecosystem_codes == code
            n_cells_ecosystem = int(mask.sum())
            if n_cells_ecosystem == 0:
                continue

            abundance = n_cells_ecosystem / n_cells_section

            ecosystem_domain_ids = domain_id[mask]
            unique_domain_ids, first_index = np.unique(ecosystem_domain_ids, return_index=True)
            ecosystem_domain_sizes = domain_size[mask][first_index]
            n_domains = len(unique_domain_ids)
            effective_n_patches = compute_effective_n_patches(ecosystem_domain_sizes)

            n_interface_edges = compute_interface_edge_count(
                sub_graph_binary, ecosystem_codes, code
            )
            mixing_index = compute_mixing_index(sub_graph_binary, ecosystem_codes, code)

            compactness_values = []
            compactness_weights = []
            for uid, size in zip(unique_domain_ids, ecosystem_domain_sizes):
                if size < 3:
                    continue
                domain_mask = ecosystem_domain_ids == uid
                iq = compute_isoperimetric_quotient(x[mask][domain_mask], y[mask][domain_mask])
                if iq is not None:
                    compactness_values.append(iq)
                    compactness_weights.append(size)
            mean_isoperimetric_quotient = (
                float(np.average(compactness_values, weights=compactness_weights))
                if compactness_values
                else np.nan
            )

            ecosystem_tumour = section_tumour.iloc[mask]
            # Report the median alongside the mean: this field is heavily
            # right-skewed (Tumour Epithelium Characterisation's tumour mask is a tight,
            # point-buffered union that "hugs" individual cell clusters
            # rather than one smoothed macro-scale margin, so a small
            # fraction of far-outlying cells inflate the mean well above
            # the typical cell's distance -- reporting mean alone would be
            # a materially misleading summary here). Some sections have no
            # tumour mask at all (`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`, `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s already-documented
            # P01_run1 finding: too few epithelial cells for any tumour
            # region to survive size-filtering) -- signed_distance is
            # legitimately all-NaN there, so mean/median are computed only
            # when at least one value exists, rather than relying on
            # pandas' silent-but-warning-emitting all-NaN reduction.
            signed_distance = ecosystem_tumour["signed_distance_to_tumour_boundary_um"]
            if signed_distance.notna().any():
                mean_signed_distance_to_tumour_boundary_um = float(signed_distance.mean())
                median_signed_distance_to_tumour_boundary_um = float(signed_distance.median())
            else:
                mean_signed_distance_to_tumour_boundary_um = np.nan
                median_signed_distance_to_tumour_boundary_um = np.nan
            fraction_inside_tumour_region = float(
                ecosystem_tumour["is_inside_tumour_region"].mean()
            )

            section_ecosystem_rows.append(
                {
                    "section_id": section_id,
                    "ecosystem_label": label,
                    "n_cells": n_cells_ecosystem,
                    "abundance": abundance,
                    "n_domains": n_domains,
                    "effective_n_patches": effective_n_patches,
                    "n_interface_edges": n_interface_edges,
                    "mixing_index": mixing_index,
                    "mean_isoperimetric_quotient": mean_isoperimetric_quotient,
                    "n_domains_with_compactness": len(compactness_values),
                    "mean_signed_distance_to_tumour_boundary_um": mean_signed_distance_to_tumour_boundary_um,
                    "median_signed_distance_to_tumour_boundary_um": median_signed_distance_to_tumour_boundary_um,
                    "fraction_inside_tumour_region": fraction_inside_tumour_region,
                }
            )

    result = pd.DataFrame(section_ecosystem_rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_rows": len(result),
        "n_sections": result["section_id"].nunique(),
        "n_ecosystems": result["ecosystem_label"].nunique(),
        "output_path": str(output_path),
    }
