#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`

Simulates spatial point patterns with known, controllable effect sizes
and known clone-spatial relationships (real-data-grounded parameters,
importance-weighted sampling from a known kernel) for `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s
null-model calibration suite. See
src/xenium_tcr_ecology/graphs/synthetic_patterns.py's module docstring

parameter grounding.

Primary output: data/graphs/synthetic/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.synthetic_patterns import build_synthetic_ground_truth_patterns


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root
        / "results"
        / "logs"
        / "09_spatial_graph_construction_and_calibration",
        script_name="07_generate_synthetic_ground_truth_patterns",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_synthetic_ground_truth_patterns(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_patients']} synthetic patient(s) x {summary['n_effect_sizes']} effect size(s) = "
        f"{summary['n_replicates_total']} replicate(s). Domain {summary['domain_size_um']}um, "
        f"length scale {summary['length_scale_um']}um. Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
