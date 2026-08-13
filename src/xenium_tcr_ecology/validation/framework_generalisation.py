"""Applies the calibrated spatial null-model framework
(`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`) end-to-end to the independent Xenium
breast-cancer dataset (`16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py`) -- the strongest test
of the registered `q1_framework_generalisation` claim (`governance/
analysis_registry.tsv`, `governance/validation_plan.tsv`, `16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py`).

Reuses `assign_tumour_region` and `sample_clone_cells`
(`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s synthetic-ground-truth generator) and
`free_permutation_pvalue`, `degree_matched_permutation_pvalue`,
`graph_distance_permutation_pvalue`, `compute_degree_strata`,
`clopper_pearson_ci` (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s null-model calibration
suite) directly, without reimplementation; the contribution here is
applying them to a different point cloud, not new statistical logic.

Uses the Xenium breast-cancer dataset's own tissue topology rather than a
second synthetic Poisson point process, since the synthetic-patient
generator's uniform Poisson background would not test whether
calibration holds on tissue with genuine spatial heterogeneity.
`assign_tumour_region`'s existing logic (a spatially contiguous region
around a random point) is applied to the Xenium breast-cancer cell
centroids, and `sample_clone_cells`'s importance-weighted sampling
injects a known synthetic clone effect on top of that topology --
observed tissue topology, known injected ground truth, since no
independent clonal ground truth exists for this dataset.

Unlike the 10 independent synthetic patients used in `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`, this
dataset is one tissue section (167,780 cells). `N_REPLICATES = 10`
independent random subsamples (`subsample_real_points`, matched to
`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s `N_CELLS_PER_PATIENT` scale) are drawn
from the same underlying point cloud, each treated as its own replicate,
mirroring the calibration suite's n=10-replicate design as closely as a
single tissue section allows -- an adaptation, not 10 independent
patients.

Comparison target: `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s established
Type I error point estimates at effect_size=0.0 (read directly from
`reports/graphs/null_model_calibration.parquet`, not re-derived), which
were themselves already elevated relative to the nominal
0.05 with wide Clopper-Pearson CIs at n=10. The
success criterion here is CI overlap with those bounds, not an exact
point-estimate match, given the calibration suite's own uncertainty at
this replicate count.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed
from xenium_tcr_ecology.graphs.candidate_graphs import build_radius_graph
from xenium_tcr_ecology.graphs.null_models import (
    NOMINAL_ALPHA,
    N_PERMUTATIONS,
    clopper_pearson_ci,
    compute_degree_strata,
    degree_matched_permutation_pvalue,
    deterministic_seed,
    free_permutation_pvalue,
    graph_distance_permutation_pvalue,
)
from xenium_tcr_ecology.graphs.synthetic_patterns import (
    EFFECT_SIZES,
    N_CELLS_PER_PATIENT,
    N_CLONE_CELLS,
    TUMOUR_FRACTION,
    assign_tumour_region,
    sample_clone_cells,
)

N_REPLICATES = 10
RNG_SEED = get_default_seed()


def subsample_real_points(points: np.ndarray, n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Spatially contiguous subsample: a random anchor point's `n_cells`
    nearest neighbours, capped at the population size.

    A uniform random subsample across this dataset's full tissue extent
    (~7,521 x 5,475 microns, 167,780 cells) destroys local spatial
    density: a uniform random 5,000-cell draw gave a mean graph degree
    of only 0.55 at the same 30um radius that gives this project's
    primary cohort a mean degree in the single digits, with 59% of
    points graph-isolated (3,788 disconnected components) -- an
    artefact of the subsampling method, not a finding about the
    framework's generalisation. A spatially contiguous k-nearest-
    neighbour window around a random anchor point preserves local
    density (matching `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s compact synthetic-patient domain
    design) and is used instead."""
    from scipy.spatial import cKDTree

    n = min(n_cells, len(points))
    anchor_idx = rng.integers(len(points))
    tree = cKDTree(points)
    _, idx = tree.query(points[anchor_idx], k=n)
    return points[idx]


def build_framework_generalisation_test(project_root: Path) -> dict:
    dataset_dir = (
        project_root / "data" / "external" / "spatial" / "Xenium_Janesick_BreastCancer_Rep1"
    )
    cells_path = dataset_dir / "Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet"
    output_path = project_root / "data" / "derived" / "framework_generalisation_results.parquet"
    return _run_framework_generalisation_test(
        project_root,
        cells_path,
        output_path,
        "`16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py`",
    )


def build_framework_generalisation_test_second_dataset(project_root: Path) -> dict:
    """Applies the framework-generalisation test above to a second,
    independent cancer type (colorectal, not breast, not HNSCC) --
    `16_external_validation_and_generalisation/09_validate_framework_on_second_cancer_type.py`,
    strengthening the `q1_framework_generalisation` claim by testing
    whether calibration holds across more than one independent tissue.
    Reuses every helper function above unchanged; only the input
    dataset and output path differ."""
    dataset_dir = (
        project_root / "data" / "external" / "spatial" / "Xenium_Oliveira_ColorectalCancer_P1"
    )
    cells_path = dataset_dir / "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cells.parquet"
    output_path = (
        project_root
        / "data"
        / "derived"
        / "framework_generalisation_results_second_dataset.parquet"
    )
    return _run_framework_generalisation_test(
        project_root,
        cells_path,
        output_path,
        "`16_external_validation_and_generalisation/08_acquire_second_independent_spatial_dataset.py`",
    )


def _run_framework_generalisation_test(
    project_root: Path, cells_path: Path, output_path: Path, acquisition_phase: str
) -> dict:
    graph_params_path = project_root / "config" / "graph_parameters.yaml"
    established_calibration_path = (
        project_root / "reports" / "graphs" / "null_model_calibration.parquet"
    )

    for path, phase in [
        (cells_path, acquisition_phase),
        (
            graph_params_path,
            "`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`",
        ),
        (
            established_calibration_path,
            "`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    radius_um = yaml.safe_load(graph_params_path.read_text())["calibrated_radius_um"]
    cells = pd.read_parquet(cells_path)
    all_points = cells[["x_centroid", "y_centroid"]].to_numpy()

    established = pd.read_parquet(established_calibration_path)
    established_null = established[established["effect_size"] == 0.0]
    established_bounds = {}
    for col in [
        "pvalue_constrained_permutation",
        "pvalue_degree_preserving",
        "pvalue_graph_preserving",
    ]:
        n_rejected = int((established_null[col] < NOMINAL_ALPHA).sum())
        established_bounds[col] = clopper_pearson_ci(n_rejected, len(established_null))

    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for replicate in range(N_REPLICATES):
        points = subsample_real_points(all_points, N_CELLS_PER_PATIENT, rng)
        is_tumour = assign_tumour_region(points, TUMOUR_FRACTION, rng)
        graph = build_radius_graph(points, radius_um)
        degree_strata = compute_degree_strata(points, radius_um)

        for effect_size in EFFECT_SIZES:
            is_clone = sample_clone_cells(
                points, is_tumour, N_CLONE_CELLS, effect_size, radius_um, rng
            )
            replicate_rng = np.random.default_rng(
                deterministic_seed(replicate, effect_size, "generalisation")
            )

            p_free = free_permutation_pvalue(
                points, is_tumour, is_clone, replicate_rng, n_permutations=N_PERMUTATIONS
            )
            p_degree = degree_matched_permutation_pvalue(
                points,
                is_tumour,
                is_clone,
                degree_strata,
                replicate_rng,
                n_permutations=N_PERMUTATIONS,
            )
            p_graph = graph_distance_permutation_pvalue(
                graph, is_tumour, is_clone, replicate_rng, n_permutations=N_PERMUTATIONS
            )

            rows.append(
                {
                    "replicate": replicate,
                    "effect_size": effect_size,
                    "pvalue_constrained_permutation": p_free,
                    "pvalue_degree_preserving": p_degree,
                    "pvalue_graph_preserving": p_graph,
                }
            )

    results = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(output_path)

    null_model_cols = [
        "pvalue_constrained_permutation",
        "pvalue_degree_preserving",
        "pvalue_graph_preserving",
    ]
    calibration_summary = {}
    for col in null_model_cols:
        by_effect = {}
        for effect_size, group in results.groupby("effect_size"):
            n_rejected = int((group[col] < NOMINAL_ALPHA).sum())
            n_total = len(group)
            ci = clopper_pearson_ci(n_rejected, n_total)
            by_effect[str(effect_size)] = {
                "rejection_rate": round(n_rejected / n_total, 4),
                "n_rejected": n_rejected,
                "n_total": n_total,
                "ci_95": [round(ci[0], 4), round(ci[1], 4)],
            }
        calibration_summary[col] = by_effect

    ci_overlap = {}
    for col in null_model_cols:
        new_ci = calibration_summary[col]["0.0"]["ci_95"]
        est_ci = established_bounds[col]
        overlaps = new_ci[0] <= est_ci[1] and est_ci[0] <= new_ci[1]
        ci_overlap[col] = overlaps

    return {
        "n_replicates": N_REPLICATES,
        "n_permutations": N_PERMUTATIONS,
        "calibration_summary": calibration_summary,
        "established_bounds_ci95": {k: list(v) for k, v in established_bounds.items()},
        "ci_overlap_with_established": ci_overlap,
        "n_null_models_overlapping": int(sum(ci_overlap.values())),
        "output_path": str(output_path),
    }
