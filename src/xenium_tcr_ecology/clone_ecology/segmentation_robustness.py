"""Segmentation robustness check (`13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py`).

The registered `segmentation_robustness_check` (governance/analysis_
registry.tsv) asks whether confirmatory clone-tumour engagement results
are robust to an independently resegmented reference subset (Phase
4.05's 3 representative sections: `P19_run1`, `P12_run2`, `P09_run2`,
morphology-free nearest-nucleus-centroid transcript reassignment).

**Scope decision:** a full recomputation of clone-tumour engagement
under the alternative resegmentation would require re-deriving cell-type
annotations (Cell Type Annotation), TCR probe calls and clone assignments (TCR Clonal Analysis),
and a spatial graph (Spatial Graph Construction and Calibration) from the resegmented expression
matrix -- effectively re-running three major phases' worth of calibrated,
multi-step pipeline on alternative data. `04_quality_control/05_resegment_reference_subset.py` itself did not
attempt this (its scope was transcript-reassignment concordance only,
`reports/qc/resegmentation_concordance.tsv`); this milestone does not
either, rather than force a rushed, under-validated partial re-
derivation that would not be a scientifically meaningful comparison
anyway. The check implemented here instead: (1) whether the primary
segmentation's engagement estimates for these 3 sections look typical
of the broader cohort (not an unusual/outlier subset within the primary
analysis itself), and (2) `04_quality_control/05_resegment_reference_subset.py`'s transcript-level concordance for these
same 3 sections, reported as an upper-bound proxy on how much a
resegmented estimate could plausibly differ, not a substitute for a
direct estimate comparison.

**Q3 barrier-topology component, closed out after `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R` (updated
here, not left permanently deferred):** the same proxy check -- (1)
above, applied to `11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s per-clone
`fibroblast_barrier_fraction`/`suppressive_myeloid_barrier_fraction`
(the two covariates `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s prespecified
`q3_barrier_topology_confirmatory` model actually uses) -- is reused
directly via `compare_resegmented_vs_rest` (the same generic,
already-tested function, no new logic needed) once `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R` existed to define
which barrier metrics are the confirmatory ones.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

RESEGMENTED_SECTIONS = ["P19_run1", "P12_run2", "P09_run2"]
ENGAGEMENT_METRICS = ["malignant_adjacency", "engagement_ratio", "fraction_inside_tumour"]
BARRIER_METRICS = ["fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction"]


def compare_resegmented_vs_rest(
    engagement: pd.DataFrame,
    resegmented_sections: list[str] = RESEGMENTED_SECTIONS,
    metrics: list[str] = ENGAGEMENT_METRICS,
) -> pd.DataFrame:
    """Pure, testable: per-metric mean comparison between clone-section
    rows inside vs. outside the resegmented reference subset."""
    in_subset = engagement[engagement["section_id"].isin(resegmented_sections)]
    outside_subset = engagement[~engagement["section_id"].isin(resegmented_sections)]
    rows = []
    for metric in metrics:
        rows.append(
            {
                "metric": metric,
                "n_in_subset": int(in_subset[metric].notna().sum()),
                "mean_in_subset": float(in_subset[metric].mean()),
                "n_outside_subset": int(outside_subset[metric].notna().sum()),
                "mean_outside_subset": float(outside_subset[metric].mean()),
            }
        )
    return pd.DataFrame(rows)


def build_segmentation_robustness(project_root: Path) -> dict:
    engagement_path = project_root / "data" / "derived" / "clone_tumour_engagement.parquet"
    barrier_path = project_root / "data" / "derived" / "clone_barrier_metrics.parquet"
    concordance_path = project_root / "reports" / "qc" / "resegmentation_concordance.tsv"
    output_path = project_root / "data" / "derived" / "segmentation_robustness_results.parquet"
    barrier_output_path = (
        project_root / "data" / "derived" / "segmentation_robustness_barrier_results.parquet"
    )

    for path, phase in [
        (engagement_path, "`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`"),
        (
            barrier_path,
            "`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`",
        ),
        (concordance_path, "`04_quality_control/05_resegment_reference_subset.py`"),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    engagement = pd.read_parquet(engagement_path)
    barrier = pd.read_parquet(barrier_path)
    concordance = pd.read_csv(concordance_path, sep="\t")

    missing_sections = set(RESEGMENTED_SECTIONS) - set(concordance["section_id"])
    if missing_sections:
        raise PipelineError(
            f"Resegmented section(s) {missing_sections} not found in '{concordance_path}'."
        )

    comparison = compare_resegmented_vs_rest(engagement, RESEGMENTED_SECTIONS, ENGAGEMENT_METRICS)
    comparison["abs_difference"] = (
        comparison["mean_in_subset"] - comparison["mean_outside_subset"]
    ).abs()
    comparison["same_direction_as_expected"] = np.sign(comparison["mean_in_subset"]) == np.sign(
        comparison["mean_outside_subset"]
    )

    barrier_comparison = compare_resegmented_vs_rest(barrier, RESEGMENTED_SECTIONS, BARRIER_METRICS)
    barrier_comparison["abs_difference"] = (
        barrier_comparison["mean_in_subset"] - barrier_comparison["mean_outside_subset"]
    ).abs()
    barrier_comparison["same_direction_as_expected"] = np.sign(
        barrier_comparison["mean_in_subset"]
    ) == np.sign(barrier_comparison["mean_outside_subset"])

    concordance_subset = concordance[concordance["section_id"].isin(RESEGMENTED_SECTIONS)][
        ["section_id", "fraction_concordant_same_cell", "pseudobulk_correlation"]
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_parquet(output_path)
    barrier_comparison.to_parquet(barrier_output_path)
    concordance_output_path = (
        project_root / "data" / "derived" / "segmentation_robustness_concordance.parquet"
    )
    concordance_subset.to_parquet(concordance_output_path)

    return {
        "n_resegmented_sections": len(RESEGMENTED_SECTIONS),
        "mean_transcript_concordance": float(
            concordance_subset["fraction_concordant_same_cell"].mean()
        ),
        "n_metrics_same_direction": int(comparison["same_direction_as_expected"].sum()),
        "n_metrics_total": len(comparison),
        "n_barrier_metrics_same_direction": int(
            barrier_comparison["same_direction_as_expected"].sum()
        ),
        "n_barrier_metrics_total": len(barrier_comparison),
        "output_path": str(output_path),
        "barrier_output_path": str(barrier_output_path),
        "concordance_output_path": str(concordance_output_path),
    }
