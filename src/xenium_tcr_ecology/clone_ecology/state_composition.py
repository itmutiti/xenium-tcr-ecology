"""Clone cell-state composition with uncertainty (`11_clone_spatial_descriptors/01_compute_clone_cell_state_composition.py`).

Per (clone, section) -- the same unit of analysis established in
`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py` (`governance/analysis_registry.tsv`'s prespecified "clone,
nested in patient and section") -- summarises the fraction of a clone's
cells in each of `06_cell_type_annotation/04_resolve_t_cell_substates.R`'s T-cell states, each with a
per-state Clopper-Pearson exact binomial confidence interval (reusing
`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s `clopper_pearson_ci`, `xenium_tcr_ecology.graphs.
null_models`, not reimplemented) rather than a bare point estimate --
most clone-section groups are small samples (`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`'s finding:
56.7% have <5 cells), where a point fraction alone would overstate
precision.

`06_cell_type_annotation/04_resolve_t_cell_substates.R`'s state taxonomy is `{Cytotoxic, Exhausted, Cycling, Treg,
CD4, CD8, Ambiguous}` -- reused exactly as defined, not remapped to the
scaffold's originally-anticipated five-category wording ("cytotoxic,
exhausted, activated, proliferative, ambiguous"): `Cycling` is the
proliferative state, there is no separate `Activated` state in this
project's taxonomy, and `Treg`/`CD4`/`CD8` are additional states beyond
those four. The taxonomy is reported as built, not remapped to force a
five-category match.

Restricted to the primary HNSCC cohort and TCR Clonal Analysis's high-confidence
clone release, identically to `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`:
every one of the 7,047 primary-cohort, release clonal cells
has a `06_cell_type_annotation/04_resolve_t_cell_substates.R` `t_cell_state`, so no cell is silently dropped for missing
state data here.

A single summary of overall composition uncertainty is also reported
per clone-section: Shannon entropy of the state composition (bits),
complementing the per-state CIs -- low entropy means a clone-section is
phenotypically homogeneous (concentrated in one state); high entropy
means it is phenotypically mixed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.clone_ecology.spatial_descriptors import filter_to_primary_cohort
from xenium_tcr_ecology.graphs.null_models import clopper_pearson_ci
from xenium_tcr_ecology.infra.exceptions import PipelineError

STATE_ORDER = ["Cytotoxic", "Exhausted", "Cycling", "Treg", "CD4", "CD8", "Ambiguous"]


def compute_state_composition_with_ci(
    states: pd.Series, state_order: list[str] = STATE_ORDER, alpha: float = 0.05
) -> dict:
    """Pure, testable: per-state fraction + Clopper-Pearson CI for one
    clone-section's cell states."""
    n = len(states)
    result = {}
    for state in state_order:
        n_successes = int((states == state).sum())
        fraction = n_successes / n if n > 0 else float("nan")
        lower, upper = clopper_pearson_ci(n_successes, n)
        result[state] = {"fraction": fraction, "ci_low": lower, "ci_high": upper}
    return result


def compute_shannon_entropy(fractions: list[float]) -> float:
    """Pure, testable: Shannon entropy (bits) of a composition vector.
    Zero-probability states contribute 0, not NaN/-inf."""
    values = np.asarray(fractions, dtype=float)
    values = values[values > 0]
    if len(values) == 0:
        return float("nan")
    return float(-np.sum(values * np.log2(values)))


def build_clone_state_composition(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    t_cell_states_path = project_root / "data" / "derived" / "t_cell_states.parquet"
    output_path = project_root / "data" / "derived" / "clone_state_composition.parquet"

    for path, phase in [
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
        (sample_manifest_path, None),
        (t_cell_states_path, "`06_cell_type_annotation/04_resolve_t_cell_substates.R`"),
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
            "No primary-cohort, release clone cells found -- cannot build state composition."
        )

    t_cell_states = pd.read_parquet(t_cell_states_path).set_index("cell_id")
    clonal["t_cell_state"] = t_cell_states.reindex(clonal.index)["t_cell_state"]
    n_missing_state = int(clonal["t_cell_state"].isna().sum())
    clonal = clonal[clonal["t_cell_state"].notna()]

    rows = []
    for (clone_id, section_id), group in clonal.groupby(["clone_id", "section_id"], observed=True):
        composition = compute_state_composition_with_ci(group["t_cell_state"])
        fractions = [composition[s]["fraction"] for s in STATE_ORDER]
        dominant_state = STATE_ORDER[int(np.argmax(fractions))] if len(group) > 0 else None

        row = {
            "clone_id": clone_id,
            "section_id": section_id,
            "patient_id": group["patient_id"].iloc[0],
            "n_cells": len(group),
            "dominant_state": dominant_state,
            "shannon_entropy_bits": compute_shannon_entropy(fractions),
        }
        for state in STATE_ORDER:
            row[f"{state.lower()}_fraction"] = composition[state]["fraction"]
            row[f"{state.lower()}_ci_low"] = composition[state]["ci_low"]
            row[f"{state.lower()}_ci_high"] = composition[state]["ci_high"]
        rows.append(row)

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_clone_section_rows": len(result),
        "n_distinct_clones": int(result["clone_id"].nunique()),
        "n_missing_state_cells_excluded": n_missing_state,
        "mean_shannon_entropy_bits": float(result["shannon_entropy_bits"].mean()),
        "output_path": str(output_path),
    }
