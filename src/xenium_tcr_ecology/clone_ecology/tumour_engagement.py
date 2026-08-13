"""Clone-tumour engagement (`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`).

The central clone-tumour-engagement measure this project exists to
characterise (following on from the `Epithelial_Tumour`/`T_cell`
local-exclusion finding in Niche and Ecosystem Discovery). Per (clone, section) -- `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`'s established unit --
four distinct quantities:

1. **Malignant-cell adjacency**: reused directly from `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`'s
   `radius_30.0um__Epithelial_Tumour` composition column (the calibrated
   30um scale, `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`) -- not recomputed via a new graph traversal,
   since `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py` already computed exactly this quantity (fraction of
   a cell's spatial neighbours that are `Epithelial_Tumour`) for every
   cell in the dataset. A clone's mean of this value across its cells is
   its malignant-cell adjacency.
2. **Contact normalised by opportunity**: raw adjacency alone conflates
   engagement with simple local tumour abundance (a clone deep in pure
   stroma has near-zero adjacency simply because no tumour is nearby,
   not because it is "avoiding" tumour) -- normalised as an enrichment
   ratio = clone's mean adjacency / that section's mean
   `Epithelial_Tumour` neighbour fraction across all its cells (the
   "opportunity" baseline), the same enrichment-ratio convention already
   used in `10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py` and `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`'s border enrichment.
3. **Penetration**: fraction of a clone's cells with negative
   `signed_distance_to_tumour_boundary_um` (`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py` -- literally
   inside the tumour mask) and the clone's minimum (deepest) signed
   distance.
4. **Interface localisation**: `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`'s `border_enrichment_
   ratio` column, joined directly by (clone_id, section_id) -- not
   redefined a second time under a different name.

**Cross-checked against Quality Control's segmentation-artefact controls, not
assumed clean:**
- `04_quality_control/06_define_qc_thresholds_hierarchically.R`'s `spillover_risk.parquet`: mean `spillover_risk_score` and
  fraction of a clone's cells flagged
  `is_boundary_adjacent_to_different_type` -- a caveat column, not a
  correction/adjustment (the flag records any different-type boundary
  adjacency, not specifically tumour-adjacency, so it cannot be
  subtracted out; it is reported alongside the engagement measure
  itself, as a segmentation-artefact-risk caveat).
- `04_quality_control/05_resegment_reference_subset.py`'s `reports/qc/resegmentation_concordance.tsv`: transcript-
  reassignment concordance under an independent, morphology-free
  resegmentation, but only available for its own 3 representative
  sections (`P19_run1`, `P12_run2`, `P09_run2`) -- reported as `NaN`,
  not fabricated, for the other 14 primary-cohort sections with no
  resegmentation check run.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.clone_ecology.spatial_descriptors import filter_to_primary_cohort
from xenium_tcr_ecology.infra.exceptions import PipelineError

TUMOUR_ADJACENCY_COLUMN = "radius_30.0um__Epithelial_Tumour"


def compute_engagement_ratio(clone_mean_adjacency: float, section_mean_adjacency: float) -> float:
    """Pure, testable: clone's malignant-cell adjacency normalised by
    the section's opportunity baseline."""
    if (
        not section_mean_adjacency
        or np.isnan(section_mean_adjacency)
        or section_mean_adjacency == 0
    ):
        return float("nan")
    return float(clone_mean_adjacency / section_mean_adjacency)


def compute_penetration(signed_distances: pd.Series) -> dict:
    """Pure, testable: penetration summary from a clone's signed
    distances to the tumour boundary (negative = inside)."""
    valid = signed_distances.dropna()
    if len(valid) == 0:
        return {"fraction_inside_tumour": float("nan"), "min_signed_distance_um": float("nan")}
    return {
        "fraction_inside_tumour": float((valid < 0).mean()),
        "min_signed_distance_um": float(valid.min()),
    }


def build_clone_tumour_engagement(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    local_compositions_path = project_root / "data" / "derived" / "local_compositions.parquet"
    tumour_boundaries_path = project_root / "data" / "derived" / "tumour_boundaries.parquet"
    spatial_descriptors_path = (
        project_root / "data" / "derived" / "clone_spatial_descriptors.parquet"
    )
    spillover_risk_path = project_root / "data" / "derived" / "spillover_risk.parquet"
    resegmentation_report_path = project_root / "reports" / "qc" / "resegmentation_concordance.tsv"
    output_path = project_root / "data" / "derived" / "clone_tumour_engagement.parquet"

    for path, phase in [
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
        (sample_manifest_path, None),
        (
            local_compositions_path,
            "`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`",
        ),
        (
            tumour_boundaries_path,
            "`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`",
        ),
        (
            spatial_descriptors_path,
            "`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`",
        ),
        (spillover_risk_path, "`04_quality_control/06_define_qc_thresholds_hierarchically.R`"),
        (resegmentation_report_path, "`04_quality_control/05_resegment_reference_subset.py`"),
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
            "No primary-cohort, release clone cells found -- cannot build tumour engagement."
        )

    local_compositions = pd.read_parquet(local_compositions_path)
    tumour_boundaries = pd.read_parquet(tumour_boundaries_path)
    spillover_risk = pd.read_parquet(spillover_risk_path)
    resegmentation = pd.read_csv(resegmentation_report_path, sep="\t").set_index("section_id")
    border_enrichment = pd.read_parquet(spatial_descriptors_path).set_index(
        ["clone_id", "section_id"]
    )["border_enrichment_ratio"]

    section_mean_adjacency = local_compositions.groupby("section_id")[
        TUMOUR_ADJACENCY_COLUMN
    ].mean()

    rows = []
    for (clone_id, section_id), group in clonal.groupby(["clone_id", "section_id"], observed=True):
        cell_ids = pd.Index(group.index)

        clone_adjacency = local_compositions.reindex(cell_ids)[TUMOUR_ADJACENCY_COLUMN]
        clone_mean_adjacency = (
            float(clone_adjacency.mean()) if clone_adjacency.notna().any() else float("nan")
        )
        engagement_ratio = compute_engagement_ratio(
            clone_mean_adjacency, section_mean_adjacency.get(section_id, float("nan"))
        )

        clone_signed_distance = tumour_boundaries.reindex(cell_ids)[
            "signed_distance_to_tumour_boundary_um"
        ]
        penetration = compute_penetration(clone_signed_distance)

        clone_spillover = spillover_risk.reindex(cell_ids)
        mean_spillover_risk_score = (
            float(clone_spillover["spillover_risk_score"].mean())
            if clone_spillover["spillover_risk_score"].notna().any()
            else float("nan")
        )
        fraction_boundary_adjacent = (
            float(clone_spillover["is_boundary_adjacent_to_different_type"].mean())
            if len(clone_spillover) > 0
            else float("nan")
        )

        has_resegmentation_check = section_id in resegmentation.index
        resegmentation_concordance = (
            float(resegmentation.loc[section_id, "fraction_concordant_same_cell"])
            if has_resegmentation_check
            else float("nan")
        )

        rows.append(
            {
                "clone_id": clone_id,
                "section_id": section_id,
                "patient_id": group["patient_id"].iloc[0],
                "n_cells": len(group),
                "malignant_adjacency": clone_mean_adjacency,
                "engagement_ratio": engagement_ratio,
                "fraction_inside_tumour": penetration["fraction_inside_tumour"],
                "min_signed_distance_um": penetration["min_signed_distance_um"],
                "interface_localisation_ratio": border_enrichment.get(
                    (clone_id, section_id), float("nan")
                ),
                "mean_spillover_risk_score": mean_spillover_risk_score,
                "fraction_boundary_adjacent_to_different_type": fraction_boundary_adjacent,
                "has_resegmentation_check": has_resegmentation_check,
                "section_resegmentation_concordance": resegmentation_concordance,
            }
        )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_clone_section_rows": len(result),
        "n_distinct_clones": int(result["clone_id"].nunique()),
        "n_rows_with_resegmentation_check": int(result["has_resegmentation_check"].sum()),
        "mean_engagement_ratio": float(result["engagement_ratio"].mean()),
        "n_clone_sections_penetrating_tumour": int((result["fraction_inside_tumour"] > 0).sum()),
        "output_path": str(output_path),
    }
