#!/usr/bin/env python3
"""
`17_statistical_closure_and_release/11_build_redesigned_manuscript_figures.py`

Builds four composite manuscript figures from already-computed,
already-validated data (no new analysis is performed here), covering
findings added after the original six main figures were finalised (the
second framework-generalisation dataset, the Q2/Q3 sensitivity checks,
and the TCR probe-vs-VDJ validation, which had no figure). See
`src/xenium_tcr_ecology/release/redesigned_main_figures.py` for the full
redesigned figure sequence and `docs/analysis_amendments.md` for why it
changed. `17_statistical_closure_and_release/02_generate_all_main_
figures.py` assembles these alongside the unchanged Figure 3 and the new
TCR validation figure.

Primary output: reports/manuscript_figures/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.release.redesigned_main_figures import (
    build_figure_1_framework_generalisation,
    build_figure_barrier_topology_with_ablation,
    build_figure_hpv_consolidated,
    build_figure_variance_partition_with_sensitivity,
)


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
        script_name="11_build_redesigned_manuscript_figures",
        project_root=project_root,
        phase="17_statistical_closure_and_release",
    )

    try:
        results = {
            "framework_generalisation": build_figure_1_framework_generalisation(project_root),
            "variance_partition_with_sensitivity": build_figure_variance_partition_with_sensitivity(
                project_root
            ),
            "barrier_topology_with_ablation": build_figure_barrier_topology_with_ablation(
                project_root
            ),
            "hpv_consolidated": build_figure_hpv_consolidated(project_root),
        }
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**{k: v["output_path"] for k, v in results.items()})
    logger.write(status="ok")
    for name, result in results.items():
        print(f"[OK]   {name}: wrote {result['output_path']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
