#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/06_compare_with_source_paper_results.py`

Compares this project's completed findings against every one of the 11
already-catalogued McCord et al. source-paper claims, and states which
findings are new. See
src/xenium_tcr_ecology/validation/source_paper_comparison.py's module
docstring.

Primary output: reports/validation/source_paper_comparison.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.source_paper_comparison import build_source_paper_comparison


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
        script_name="06_compare_with_source_paper_results",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_source_paper_comparison(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_claims_compared']} source-paper claim(s) compared: {summary['status_counts']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
