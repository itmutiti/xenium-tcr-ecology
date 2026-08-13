#!/usr/bin/env python3
"""
`01_project_setup_and_governance/00_validate_project_scope.py`

Creates project_charter.yaml, separating published findings (McCord et al.
2026) from this project's original hypotheses (Q1-Q3).

Always succeeds -- this is a factual record, not an approval gate.

Primary output: project_charter.yaml
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.governance.charter import build_project_charter
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    output_path = project_root / "project_charter.yaml"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "01_project_setup_and_governance",
        script_name="00_validate_project_scope",
        project_root=project_root,
        phase="01_project_setup_and_governance",
    )

    charter = build_project_charter(project_root, output_path)

    logger.log_event(
        primary_questions=charter["primary_questions"],
        output=str(output_path),
    )
    logger.write(status="ok")
    print(f"[OK]   Wrote {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
