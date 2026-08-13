"""Synthetic ground-truth spatial pattern generation (`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`).

Simulates spatial point patterns with known, controllable effect sizes and
known clone-spatial relationships, for `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s null-model calibration
suite (the registered `null_model_calibration_suite` analysis,
`governance/analysis_registry.tsv`) to test Type I error and power
against, before any null model is trusted on real clones.

**Parameters, grounded directly in this project's own data, not
arbitrary:** `N_SYNTHETIC_PATIENTS = 10` (the registered analysis'
explicit "n=10-patient scale" target). `TUMOUR_FRACTION = 0.073` (the
fraction of all cells inside a validated tumour region across this
project's actual 18 sections, `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s `tumour_boundaries.parquet`
-- not an arbitrary round number). `N_CLONE_CELLS = 30` (near this
project's 75th-percentile clone size, `08_tcr_clonal_analysis/07_build_clone_metadata_table.py`'s
`clone_metadata.parquet`, giving the calibration suite a realistic,
not artificially generous, amount of signal to detect). `DOMAIN_SIZE_UM`
derived from `N_CELLS_PER_PATIENT` and the mean cell area (Phase
4.00, ~90 square-um/cell) at a plausible ~70% tissue packing fraction.

**Generative model, a standard, well-understood method (importance-
weighted sampling from a known kernel), not an ad hoc heuristic:**
1. Background cells: a uniform Poisson point process in a square domain.
2. Tumour region: a spatially contiguous subset (all background points
   within a radius of one random centre, sized to hit
   `TUMOUR_FRACTION`) -- not randomly scattered tumour labels, which
   would not resemble tissue architecture (tumour masses are
   contiguous, `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`, `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`).
3. Clone cells: sampled without replacement from the remaining
   (non-tumour) points, with each point's selection weight
   `exp(-effect_size * distance_to_nearest_tumour_cell / length_scale)`
   -- `effect_size = 0.0` gives uniform (i.e. unweighted, uniform random)
   sampling, a true null with no spatial association by construction;
   `effect_size > 0` gives a controllable, known attraction toward the
   tumour region, strictly increasing in strength. `length_scale` is
   `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`'s calibrated radius (30um), so the effect size is
   expressed in units directly comparable to this project's calibrated
   spatial scale.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree

from xenium_tcr_ecology.infra.exceptions import PipelineError

N_SYNTHETIC_PATIENTS = 10
N_CELLS_PER_PATIENT = 5000
TUMOUR_FRACTION = 0.073
N_CLONE_CELLS = 30
EFFECT_SIZES = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 1.0]
MEAN_CELL_AREA_UM2 = 90.0
TISSUE_PACKING_FRACTION = 0.7


def compute_domain_size_um(n_cells: int = N_CELLS_PER_PATIENT) -> float:
    total_area = n_cells * MEAN_CELL_AREA_UM2 / TISSUE_PACKING_FRACTION
    return float(np.sqrt(total_area))


def generate_background_points(
    n_cells: int, domain_size_um: float, rng: np.random.Generator
) -> np.ndarray:
    return rng.uniform(0, domain_size_um, size=(n_cells, 2))


def assign_tumour_region(
    points: np.ndarray, tumour_fraction: float, rng: np.random.Generator
) -> np.ndarray:
    """Pure, testable contiguous tumour-region assignment: all points
    within a radius of one random centre, sized to hit `tumour_fraction`
    of the population."""
    n = len(points)
    n_tumour = max(1, int(round(n * tumour_fraction)))
    centre = points[rng.integers(n)]
    distances = np.linalg.norm(points - centre, axis=1)
    tumour_indices = np.argsort(distances)[:n_tumour]
    mask = np.zeros(n, dtype=bool)
    mask[tumour_indices] = True
    return mask


def sample_clone_cells(
    points: np.ndarray,
    is_tumour: np.ndarray,
    n_clone_cells: int,
    effect_size: float,
    length_scale: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Pure, testable importance-weighted clone-cell sampling. Returns a
    boolean mask the same length as `points`. `effect_size = 0.0` gives
    uniform random sampling (a true null); larger values give
    increasing, known attraction toward the tumour region."""
    n = len(points)
    candidate_mask = ~is_tumour
    candidate_indices = np.where(candidate_mask)[0]
    n_clone_cells = min(n_clone_cells, len(candidate_indices))

    if effect_size == 0.0 or not is_tumour.any():
        chosen = rng.choice(candidate_indices, size=n_clone_cells, replace=False)
    else:
        tumour_tree = cKDTree(points[is_tumour])
        distances, _ = tumour_tree.query(points[candidate_indices])
        weights = np.exp(-effect_size * distances / length_scale)
        weights = weights / weights.sum()
        chosen = rng.choice(candidate_indices, size=n_clone_cells, replace=False, p=weights)

    mask = np.zeros(n, dtype=bool)
    mask[chosen] = True
    return mask


def build_synthetic_ground_truth_patterns(project_root: Path) -> dict:
    graph_params_path = project_root / "config" / "graph_parameters.yaml"
    output_dir = project_root / "data" / "graphs" / "synthetic"

    if not graph_params_path.is_file():
        raise PipelineError(
            f"'{graph_params_path}' not found. Run `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` first."
        )
    length_scale = yaml.safe_load(graph_params_path.read_text())["calibrated_radius_um"]

    domain_size_um = compute_domain_size_um()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for patient_idx in range(N_SYNTHETIC_PATIENTS):
        patient_id = f"SYN{patient_idx:02d}"
        for effect_size in EFFECT_SIZES:
            rng = np.random.default_rng(hash((patient_idx, effect_size)) % (2**32))
            points = generate_background_points(N_CELLS_PER_PATIENT, domain_size_um, rng)
            is_tumour = assign_tumour_region(points, TUMOUR_FRACTION, rng)
            is_clone = sample_clone_cells(
                points, is_tumour, N_CLONE_CELLS, effect_size, length_scale, rng
            )

            replicate = pd.DataFrame(
                {"x": points[:, 0], "y": points[:, 1], "is_tumour": is_tumour, "is_clone": is_clone}
            )
            safe_effect = str(effect_size).replace(".", "p")
            out_path = output_dir / f"{patient_id}_effect{safe_effect}.parquet"
            replicate.to_parquet(out_path)

            rows.append(
                {
                    "patient_id": patient_id,
                    "effect_size": effect_size,
                    "n_cells": len(replicate),
                    "n_tumour_cells": int(is_tumour.sum()),
                    "n_clone_cells": int(is_clone.sum()),
                    "output_path": str(out_path),
                }
            )

    manifest = pd.DataFrame(rows)
    manifest_path = output_dir / "_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_patients": N_SYNTHETIC_PATIENTS,
        "n_effect_sizes": len(EFFECT_SIZES),
        "n_replicates_total": len(manifest),
        "domain_size_um": round(domain_size_um, 1),
        "length_scale_um": length_scale,
        "output_dir": str(output_dir),
        "manifest_path": str(manifest_path),
    }
