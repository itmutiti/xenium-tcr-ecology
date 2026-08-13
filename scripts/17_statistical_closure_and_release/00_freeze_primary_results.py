#!/usr/bin/env python3
"""
`17_statistical_closure_and_release/00_freeze_primary_results.py`

Freezes the prespecified primary results (Q1-Q3, plus the single HPV
contrast set) into an immutable, hash-manifested release directory
before any exploratory extension work builds on them. See
src/xenium_tcr_ecology/release/freeze_primary_results.py's module
docstring.

Primary output: data/releases/final_primary/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.release.freeze_primary_results import build_primary_results_freeze


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
        script_name="00_freeze_primary_results",
        project_root=project_root,
        phase="17_statistical_closure_and_release",
    )

    try:
        summary = build_primary_results_freeze(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_primary_analyses']} prespecified primary result(s) frozen ({summary['n_files']} files). {summary['release_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
