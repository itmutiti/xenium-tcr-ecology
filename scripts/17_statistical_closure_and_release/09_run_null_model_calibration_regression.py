#!/usr/bin/env python3
"""
`17_statistical_closure_and_release/09_run_null_model_calibration_regression.py`

Re-runs `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s null-model calibration suite against the
frozen release code, guarding against calibration drift. See
src/xenium_tcr_ecology/release/calibration_regression.py's module
docstring.

Primary output: reports/release/calibration_regression.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.release.calibration_regression import build_calibration_regression


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "17_statistical_closure_and_release",
        script_name="09_run_null_model_calibration_regression",
        project_root=project_root,
        phase="17_statistical_closure_and_release",
    )

    try:
        summary = build_calibration_regression(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_ci_overlap']}/{summary['n_comparisons']} comparison(s) CI-overlap (no drift flagged). Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
