#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/09_validate_framework_on_second_cancer_type.py`

Applies this project's already-calibrated spatial null-model framework
(`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`) end-to-end to a second,
independent Xenium dataset (`08_acquire_second_independent_
spatial_dataset.py`, colorectal cancer, de Oliveira et al. 2025) --
strengthening the registered `q1_framework_generalisation` claim by
testing calibration on more than one independent tissue type, not
re-testing the breast-cancer dataset
(`05_validate_framework_on_independent_dataset.py`) a second time. See
src/xenium_tcr_ecology/validation/framework_generalisation.py's module
docstring.

Primary output: reports/validation/framework_generalisation_second_dataset.pdf
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.framework_generalisation import (
    build_framework_generalisation_test_second_dataset,
)

MODEL_LABELS = {
    "pvalue_constrained_permutation": "Constrained permutation",
    "pvalue_degree_preserving": "Degree-preserving",
    "pvalue_graph_preserving": "Graph-preserving",
}


def render_calibration_plot(calibration_summary: dict, output_path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    for col, by_effect in calibration_summary.items():
        sorted_keys = sorted(by_effect, key=float)
        effect_sizes = [float(k) for k in sorted_keys]
        rates = [by_effect[k]["rejection_rate"] for k in sorted_keys]
        ax.plot(effect_sizes, rates, marker="o", label=MODEL_LABELS.get(col, col))
    ax.axhline(0.05, linestyle="--", color="grey", label="nominal alpha=0.05")
    ax.set_xlabel("Effect size (0 = true null)")
    ax.set_ylabel("Rejection rate (p < 0.05) across 10 subsample replicates")
    ax.set_title("Framework calibration on the independent Xenium colorectal-cancer dataset")
    ax.legend()
    ax.set_ylim(-0.05, 1.05)
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "16_external_validation_and_generalisation",
        script_name="09_validate_framework_on_second_cancer_type",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_framework_generalisation_test_second_dataset(project_root)
        report_path = (
            project_root / "reports" / "validation" / "framework_generalisation_second_dataset.pdf"
        )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        render_calibration_plot(summary["calibration_summary"], report_path)
        summary["report_path"] = str(report_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_replicates']} replicate(s), {summary['n_permutations']} permutations each."
    )
    for model in summary["calibration_summary"]:
        established = summary["established_bounds_ci95"][model]
        new_ci = summary["calibration_summary"][model]["0.0"]["ci_95"]
        overlaps = summary["ci_overlap_with_established"][model]
        print(
            f"       {model}: Type I error CI={new_ci} vs. established={established} -> overlap={overlaps}"
        )
    print(
        f"[OK]   {summary['n_null_models_overlapping']}/3 null model(s) CI-overlap with `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s established bounds."
    )
    print(f"[OK]   Wrote {summary['output_path']}, {summary['report_path']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
