#!/usr/bin/env python3
"""
`17_statistical_closure_and_release/06_build_public_data_release.py`

Exports non-sensitive derived data, schemas and lightweight examples with
a licence, into a hash-manifested public-release bundle. See
`src/xenium_tcr_ecology/release/public_data_release.py`'s module
docstring for the exact scope (three source categories) and what is
excluded.

Primary output: release/data/

This script assembles the release bundle. It does not perform the
privacy/licensing review the release procedure requires.
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.release.public_data_release import build_public_data_release


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
        script_name="06_build_public_data_release",
        project_root=project_root,
        phase="17_statistical_closure_and_release",
    )

    try:
        summary = build_public_data_release(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Public data release built: {summary['n_files']} file(s), "
        f"{summary['license']}. {summary['release_dir']}\n"
        "       This bundle still requires a privacy/licensing review "
        "before upload."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
