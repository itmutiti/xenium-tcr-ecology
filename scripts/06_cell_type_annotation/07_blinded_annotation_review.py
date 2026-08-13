#!/usr/bin/env python3
"""
`06_cell_type_annotation/07_blinded_annotation_review.py`

Generates blinded spatial review panels (a stratified, ambiguous-cell-
oversampled sample) and a correctly-structured, empty adjudication log
template. Does not fabricate adjudications: judging whether a specific
cell's tissue morphology matches its algorithmic prediction requires
a human domain expert actually looking at each panel. Completing this
phase (filling in adjudication_log.tsv) is a human-in-the-loop
step -- see src/xenium_tcr_ecology/annotation/blinded_review.py's module
docstring.

Primary output: reports/annotation/adjudication_log.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.blinded_review import build_blinded_review_panels


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "06_cell_type_annotation",
        script_name="07_blinded_annotation_review",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_blinded_review_panels(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_panels_rendered']}/{summary['n_panels_requested']} panel(s) rendered "
        f"({summary['n_panels_failed']} failed), {summary['n_ambiguous_in_sample']} from ambiguous cells. "
        f"Wrote {summary['panels_dir']}, {summary['key_path']}, {summary['adjudication_log_path']}. "
        f"STATUS: {summary['adjudication_status']} -- human expert review required to complete this phase."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
