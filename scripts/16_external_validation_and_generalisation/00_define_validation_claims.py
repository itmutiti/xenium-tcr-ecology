#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/00_define_validation_claims.py`

Predeclares this project's external-validation claims, datasets
and success criteria before any External Validation and Generalisation dataset is acquired
or analysed. See src/xenium_tcr_ecology/validation/plan.py's module
docstring.

Primary output: governance/validation_plan.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.plan import build_validation_plan


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
        script_name="00_define_validation_claims",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_validation_plan(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_claims']} validation claim(s) predeclared. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
