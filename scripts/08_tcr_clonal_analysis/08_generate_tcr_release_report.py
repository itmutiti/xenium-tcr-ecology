#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`

Freezes high-confidence clone definitions (`08_tcr_clonal_analysis/07_build_clone_metadata_table.py`'s clones with at
least one `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py` `singlet`-resolution cell) for primary analysis and
documents excluded/ambiguous calls, following the same SHA256-manifested
freeze pattern as `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`. See
src/xenium_tcr_ecology/tcr/release_freeze.py's module docstring for the high-confidence criterion and
CONDITIONAL GO release status rationale.

Primary output: data/releases/v1_tcr_calls/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.release_freeze import RELEASE_NAME, freeze_tcr_calls


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "08_tcr_clonal_analysis",
        script_name="08_generate_tcr_release_report",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    release_dir = project_root / "data" / "releases" / RELEASE_NAME

    try:
        summary = freeze_tcr_calls(project_root, release_dir)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   STATUS: {summary['status']}. {summary['n_high_confidence_clones']} high-confidence clone(s) "
        f"({summary['n_high_confidence_cells']} cells), {summary['n_excluded_clones']} excluded clone(s). "
        f"Wrote {summary['release_dir']} ({summary['n_files']} files)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
