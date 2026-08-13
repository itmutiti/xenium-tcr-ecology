#!/usr/bin/env python3
"""
`17_statistical_closure_and_release/05_render_reproducible_manuscript.qmd` -- 05b_render_manuscript_pdf.py

Companion driver for `05_render_reproducible_manuscript.qmd`: renders
the `.qmd` source content to PDF via a pure-Python `markdown` +
`xhtml2pdf` path. See
src/xenium_tcr_ecology/release/manuscript.py's module docstring.

Primary output: manuscript/manuscript.pdf
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.release.manuscript import build_manuscript_pdf


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
        script_name="05b_render_manuscript_pdf",
        project_root=project_root,
        phase="17_statistical_closure_and_release",
    )

    try:
        summary = build_manuscript_pdf(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Manuscript rendered ({summary['n_body_characters']} body characters). Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
