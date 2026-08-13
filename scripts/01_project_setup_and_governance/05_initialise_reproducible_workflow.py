#!/usr/bin/env python3
"""
`01_project_setup_and_governance/05_initialise_reproducible_workflow.py`

Validates the Snakefile, config/config.yaml, and config/global_seed.yaml
(created by the repository scaffold) and adds the local Snakemake execution
profile (profiles/local/config.yaml) if not already present. Idempotent:
re-running without --force leaves an existing profile untouched.

Primary output: Snakefile; config/config.yaml; profiles/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.governance.workflow_init import initialise_workflow
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root


def main() -> int:
    parser = base_parser(__doc__)
    parser.add_argument("--force", action="store_true", help="Overwrite an existing local profile.")
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "01_project_setup_and_governance",
        script_name="05_initialise_reproducible_workflow",
        project_root=project_root,
        phase="01_project_setup_and_governance",
    )

    try:
        summary = initialise_workflow(project_root, force=args.force)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(f"[OK]   Snakefile/config validated. default_seed={summary['default_seed']}.")
    print(
        f"[OK]   Local profile {'written' if summary['profile_written'] else 'already present (kept)'}: {summary['profile_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
