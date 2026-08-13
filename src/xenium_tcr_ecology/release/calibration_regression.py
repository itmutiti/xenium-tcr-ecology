"""Re-runs `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s null-model calibration suite against
the frozen release code, guarding against calibration drift (Phase
17.09) -- the Dockerfile's header comment names this exact script as
the standard pre-release verification step.

**Direct consequence of `17_statistical_closure_and_release/08_run_end_to_end_reproduction_test.sh`'s finding:** this
re-run is what actually exercises the `deterministic_seed` fix
(`xenium_tcr_ecology.graphs.null_models`) added after `17_statistical_closure_and_release/08_run_end_to_end_reproduction_test.sh`'s
reproducibility-gap discovery -- before that fix, re-running this exact
calibration suite in a fresh process would draw different
random permutations each time (Python's built-in `hash()` on a
`str`-containing tuple is randomised per process), making a literal
"regression against the stored result" comparison unable to distinguish
drift from ordinary re-seeding noise. Comparison here uses
Clopper-Pearson CI overlap (the same method Phase 16.05 already used
for its cross-dataset comparison), not an exact value match, since even
a fully deterministic re-run at only 10 replicates carries expected
sampling variability of its own if any upstream data changed.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.graphs.null_models import (
    NOMINAL_ALPHA,
    build_null_model_calibration,
    clopper_pearson_ci,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError

NULL_MODEL_COLUMNS = [
    "pvalue_constrained_permutation",
    "pvalue_degree_preserving",
    "pvalue_graph_preserving",
]


def ci_overlaps(ci_a: tuple[float, float], ci_b: tuple[float, float]) -> bool:
    """Pure, testable: interval-overlap check (closed intervals,
    touching endpoints count as overlapping)."""
    return ci_a[0] <= ci_b[1] and ci_b[0] <= ci_a[1]


def build_calibration_regression(project_root: Path) -> dict:
    stored_path = project_root / "reports" / "graphs" / "null_model_calibration.parquet"
    if not stored_path.exists():
        raise PipelineError(
            f"'{stored_path}' not found. Run `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py` first."
        )

    stored = pd.read_parquet(stored_path)

    fresh_summary = build_null_model_calibration(project_root)
    fresh = pd.read_parquet(fresh_summary["output_path"])

    rows = []
    for col in NULL_MODEL_COLUMNS:
        for effect_size in sorted(stored["effect_size"].unique()):
            stored_group = stored[stored["effect_size"] == effect_size]
            fresh_group = fresh[fresh["effect_size"] == effect_size]
            stored_n_rej = int((stored_group[col] < NOMINAL_ALPHA).sum())
            fresh_n_rej = int((fresh_group[col] < NOMINAL_ALPHA).sum())
            stored_ci = clopper_pearson_ci(stored_n_rej, len(stored_group))
            fresh_ci = clopper_pearson_ci(fresh_n_rej, len(fresh_group))
            overlaps = ci_overlaps(stored_ci, fresh_ci)
            rows.append(
                {
                    "null_model": col,
                    "effect_size": effect_size,
                    "stored_rejection_rate": stored_n_rej / len(stored_group),
                    "fresh_rejection_rate": fresh_n_rej / len(fresh_group),
                    "ci_overlap": overlaps,
                }
            )

    result = pd.DataFrame(rows)
    output_dir = project_root / "reports" / "release"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "calibration_regression.tsv"
    result.to_csv(output_path, sep="\t", index=False)

    n_overlap = int(result["ci_overlap"].sum())

    return {
        "n_comparisons": len(result),
        "n_ci_overlap": n_overlap,
        "n_drift_flagged": len(result) - n_overlap,
        "output_path": str(output_path),
    }
