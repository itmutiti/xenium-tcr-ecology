#!/usr/bin/env python3
"""
`17_statistical_closure_and_release/02_generate_all_main_figures.py`

Assembles the manuscript main figures from already-generated,
already-validated PDF outputs -- no analysis is regenerated here.
See src/xenium_tcr_ecology/release/main_figures.py's module docstring
.

Primary output: figures/main/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.release.main_figures import build_main_figures


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
        script_name="02_generate_all_main_figures",
        project_root=project_root,
        phase="17_statistical_closure_and_release",
    )

    try:
        summary = build_main_figures(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(f"[OK]   {summary['n_figures']} main figure(s) assembled. {summary['figures_dir']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
