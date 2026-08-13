#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/06_validate_boundaries_against_morphology.py`

Overlays `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s tumour-boundary polygon outline on the
DAPI morphology image, at points sampled along the boundary itself, and
generates a correctly-structured, empty manual-review log template. Does
not fabricate a "manual-review agreement" figure: no pathologist tumour/
normal annotation exists in this dataset (same finding
as `07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`) -- see
src/xenium_tcr_ecology/tumour/boundary_validation.py's module docstring
.

Primary output: reports/tumour/boundary_validation.pdf
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.boundary_validation import build_boundary_validation_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "07_tumour_epithelium_characterisation",
        script_name="06_validate_boundaries_against_morphology",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_boundary_validation_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_panels_rendered']} panel(s) rendered ({summary['n_panels_failed']} failed), "
        f"{summary['n_sections_no_boundary']} section(s) with no tumour boundary to render. "
        f"Wrote {summary['output_path']}, {summary['review_log_path']}. "
        f"STATUS: {summary['review_status']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
