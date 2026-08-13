"""Clone spatial descriptors, rarefaction-normalised (`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`).

Computes five spatial descriptors per (clone, section) -- the unit of
analysis (`governance/analysis_registry.tsv`'s prespecified "clone,
nested in patient and section", confirmed: 158 primary-cohort
clones span 261 distinct clone-section pairs, 103 clones present in
both of a patient's technical-replicate runs). Restricted to
`metadata/sample_manifest.tsv`'s `included_in_primary_hnscc_cohort ==
True` sections throughout (`P10_run1`'s specimen is an
Ameloblastoma, not HNSCC, and contributes 8 clones to TCR Clonal Analysis's
release that must not silently enter Clone Spatial Descriptors+).

A clone's cells are exactly `08_tcr_clonal_analysis/07_build_clone_metadata_table.py`'s definition (`resolution in
{"singlet", "low_confidence"}` cells whose `detected_probes` equals the
clone's `clone_id`), restricted to TCR Clonal Analysis's `high_confidence_clones`
release -- not an independently re-derived clone definition.

Three of the five descriptors (dispersion, clustering, domain occupancy)
are point-pattern statistics with an n-driven bias (more cells trivially
look "more spread out" / "occupy more domains" purely from sampling more
of the tissue) -- the same statistical structure as species-accumulation
curves in classical ecology, which this project's "clone ecology"
framing deliberately invokes. These three are computed via rarefaction:
repeated fixed-size (`N_RAREFY = 5`) random subsampling without
replacement, averaged over `B_RAREFY_ITERATIONS = 100` draws; clones
with fewer than `N_RAREFY` cells get `NaN` for these three (a value that
cannot be extrapolated below a clone's true size) rather than a
fabricated estimate. The other two descriptors (tumour distance, border
enrichment) are per-cell averages, which are unbiased at any sample size
(noisier for small clones, but not systematically biased the way
pattern statistics are) -- computed directly on each clone's full cell
set, not rarefied.

**Method per descriptor:**
- Dispersion: standard distance (sqrt of mean squared distance from the
  rarefied subsample's centroid).
- Clustering: Clark-Evans nearest-neighbour index R = observed / expected
  mean nearest-neighbour distance under complete spatial randomness, at
  that section's cell density (not the clone's own -- too few points for
  a meaningful self-density estimate). R<1 = more clustered than random;
  R>1 = more dispersed. Uses the standard, uncorrected 2D Poisson-process
  expectation (`0.5/sqrt(density)`) -- a documented limitation: no
  edge-effect correction is applied, which is known to bias Clark-Evans
  R toward "more regular" near a tissue section's boundary; not
  corrected for here.
- Domain occupancy: number of distinct `10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py` tissue `domain_id`
  values among the rarefied subsample.
- Tumour distance: mean and median `signed_distance_to_tumour_
  boundary_um` (`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`) across a clone's full cell set.
- Border enrichment: fraction of a clone's cells in `07_tumour_epithelium_characterisation/07_define_invasive_front_and_compartments.py`'s
  primary (3um-band) `inner_margin`/`outer_margin` compartments, divided
  by that section's inner/outer-margin fraction -- an enrichment ratio,
  matching `10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`'s enrichment-ratio naming convention rather than a
  new ad hoc definition.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError
from scipy.spatial.distance import cdist

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed

N_RAREFY = 5
B_RAREFY_ITERATIONS = 100
RNG_SEED = get_default_seed()
MARGIN_COMPARTMENTS = {"inner_margin", "outer_margin"}


def filter_to_primary_cohort(sample_manifest: pd.DataFrame) -> pd.DataFrame:
    """Pure, testable: restricts to sections flagged
    `included_in_primary_hnscc_cohort == True`."""
    return sample_manifest[sample_manifest["included_in_primary_hnscc_cohort"]]


def compute_dispersion(x: np.ndarray, y: np.ndarray) -> float:
    """Pure, testable: standard distance (spread) from the centroid."""
    cx, cy = x.mean(), y.mean()
    return float(np.sqrt(np.mean((x - cx) ** 2 + (y - cy) ** 2)))


def compute_convex_hull_area(x: np.ndarray, y: np.ndarray) -> float | None:
    """Pure, testable: convex-hull area; None for <3 points or a
    degenerate (collinear) hull, never a fabricated value."""
    if len(x) < 3:
        return None
    try:
        hull = ConvexHull(np.column_stack([x, y]))
    except QhullError:
        return None
    return float(hull.volume)


def compute_clark_evans_index(x: np.ndarray, y: np.ndarray, reference_density: float) -> float:
    """Pure, testable: observed/expected mean nearest-neighbour distance
    under complete spatial randomness at `reference_density` (cells per
    unit area from the section, not the clone's own sparse points)."""
    n = len(x)
    if n < 2 or not reference_density or reference_density <= 0 or np.isnan(reference_density):
        return float("nan")
    points = np.column_stack([x, y])
    dist = cdist(points, points)
    np.fill_diagonal(dist, np.inf)
    observed_mean_nn = float(dist.min(axis=1).mean())
    expected_mean_nn = 0.5 / np.sqrt(reference_density)
    return observed_mean_nn / expected_mean_nn


def compute_border_enrichment(
    clone_compartments: pd.Series, section_compartments: pd.Series
) -> float:
    """Pure, testable: clone's margin-compartment fraction divided by
    the section's margin-compartment fraction."""
    section_frac = section_compartments.isin(MARGIN_COMPARTMENTS).mean()
    if not section_frac or np.isnan(section_frac):
        return float("nan")
    clone_frac = clone_compartments.isin(MARGIN_COMPARTMENTS).mean()
    return float(clone_frac / section_frac)


def compute_rarefied_dispersion(
    x: np.ndarray, y: np.ndarray, n_rarefy: int, n_iterations: int, rng: np.random.Generator
) -> float:
    """Pure, testable (given a seeded `rng`): mean dispersion across
    `n_iterations` random size-`n_rarefy` subsamples. NaN if fewer than
    `n_rarefy` cells are available."""
    n = len(x)
    if n < n_rarefy:
        return float("nan")
    values = []
    for _ in range(n_iterations):
        idx = rng.choice(n, size=n_rarefy, replace=False)
        values.append(compute_dispersion(x[idx], y[idx]))
    return float(np.mean(values))


def compute_rarefied_clark_evans(
    x: np.ndarray,
    y: np.ndarray,
    reference_density: float,
    n_rarefy: int,
    n_iterations: int,
    rng: np.random.Generator,
) -> float:
    """Pure, testable (given a seeded `rng`): mean Clark-Evans index
    across `n_iterations` random size-`n_rarefy` subsamples."""
    n = len(x)
    if n < n_rarefy:
        return float("nan")
    values = []
    for _ in range(n_iterations):
        idx = rng.choice(n, size=n_rarefy, replace=False)
        values.append(compute_clark_evans_index(x[idx], y[idx], reference_density))
    return float(np.mean(values))


def compute_rarefied_domain_richness(
    domain_ids: np.ndarray, n_rarefy: int, n_iterations: int, rng: np.random.Generator
) -> float:
    """Pure, testable (given a seeded `rng`): mean number of distinct
    domains among `n_iterations` random size-`n_rarefy` subsamples."""
    n = len(domain_ids)
    if n < n_rarefy:
        return float("nan")
    values = []
    for _ in range(n_iterations):
        idx = rng.choice(n, size=n_rarefy, replace=False)
        values.append(len(np.unique(domain_ids[idx])))
    return float(np.mean(values))


def build_clone_spatial_descriptors(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    tumour_boundaries_path = project_root / "data" / "derived" / "tumour_boundaries.parquet"
    spatial_compartments_path = project_root / "data" / "derived" / "spatial_compartments.parquet"
    tissue_domains_path = project_root / "data" / "derived" / "tissue_domains.parquet"
    output_path = project_root / "data" / "derived" / "clone_spatial_descriptors.parquet"

    for path, phase in [
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
        (sample_manifest_path, None),
        (
            primary_graphs_dir,
            "`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`",
        ),
        (
            tumour_boundaries_path,
            "`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`",
        ),
        (
            spatial_compartments_path,
            "`07_tumour_epithelium_characterisation/07_define_invasive_front_and_compartments.py`",
        ),
        (tissue_domains_path, "`10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`"),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found." + (f" Run {phase} first." if phase else ""))

    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
    primary_manifest = filter_to_primary_cohort(sample_manifest)
    primary_sections = set(primary_manifest["section_id"])

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
            "No primary-cohort, release clone cells found -- cannot build spatial descriptors."
        )

    tumour_boundaries = pd.read_parquet(tumour_boundaries_path)
    spatial_compartments = pd.read_parquet(spatial_compartments_path)
    tissue_domains = pd.read_parquet(tissue_domains_path).set_index("cell_id")

    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for section_id in sorted(clonal["section_id"].unique()):
        node_metadata = pd.read_csv(
            primary_graphs_dir / section_id / "node_metadata.tsv", sep="\t"
        ).set_index("cell_id")
        section_x = node_metadata["x_centroid"].to_numpy()
        section_y = node_metadata["y_centroid"].to_numpy()
        section_area = compute_convex_hull_area(section_x, section_y)
        section_density = len(node_metadata) / section_area if section_area else float("nan")
        section_compartments = spatial_compartments.reindex(node_metadata.index)["compartment"]

        section_clonal = clonal[clonal["section_id"] == section_id]
        for clone_id, group in section_clonal.groupby("clone_id"):
            cell_ids = pd.Index(group.index)
            coords = node_metadata.reindex(cell_ids)
            valid = coords["x_centroid"].notna().to_numpy()
            cell_ids = cell_ids[valid]
            x = coords.loc[valid, "x_centroid"].to_numpy()
            y = coords.loc[valid, "y_centroid"].to_numpy()
            n_cells = len(cell_ids)
            if n_cells == 0:
                continue

            dispersion = compute_rarefied_dispersion(x, y, N_RAREFY, B_RAREFY_ITERATIONS, rng)
            clustering = compute_rarefied_clark_evans(
                x, y, section_density, N_RAREFY, B_RAREFY_ITERATIONS, rng
            )

            clone_domains = tissue_domains.reindex(cell_ids)["domain_id"]
            valid_domain = clone_domains.notna().to_numpy()
            n_distinct_domains_observed = int(clone_domains[valid_domain].nunique())
            domain_richness = (
                compute_rarefied_domain_richness(
                    clone_domains[valid_domain].to_numpy().astype(int),
                    N_RAREFY,
                    B_RAREFY_ITERATIONS,
                    rng,
                )
                if valid_domain.sum() >= N_RAREFY
                else float("nan")
            )

            clone_tumour_distance = tumour_boundaries.reindex(cell_ids)[
                "signed_distance_to_tumour_boundary_um"
            ]
            has_tumour_distance = clone_tumour_distance.notna().any()
            mean_tumour_distance = (
                float(clone_tumour_distance.mean()) if has_tumour_distance else float("nan")
            )
            median_tumour_distance = (
                float(clone_tumour_distance.median()) if has_tumour_distance else float("nan")
            )

            clone_compartments = section_compartments.reindex(cell_ids)
            border_enrichment = compute_border_enrichment(clone_compartments, section_compartments)

            rows.append(
                {
                    "clone_id": clone_id,
                    "section_id": section_id,
                    "patient_id": group["patient_id"].iloc[0],
                    "n_cells": n_cells,
                    "dispersion_rarefied": dispersion,
                    "clark_evans_index_rarefied": clustering,
                    "domain_richness_rarefied": domain_richness,
                    "n_distinct_domains_observed": n_distinct_domains_observed,
                    "mean_tumour_distance_um": mean_tumour_distance,
                    "median_tumour_distance_um": median_tumour_distance,
                    "border_enrichment_ratio": border_enrichment,
                }
            )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_clone_section_rows": len(result),
        "n_distinct_clones": int(result["clone_id"].nunique()),
        "n_sections": int(result["section_id"].nunique()),
        "n_rarefied_dispersion": int(result["dispersion_rarefied"].notna().sum()),
        "output_path": str(output_path),
    }
