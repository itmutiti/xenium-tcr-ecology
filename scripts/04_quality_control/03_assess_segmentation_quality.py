#!/usr/bin/env python3
"""
`04_quality_control/03_assess_segmentation_quality.py`

Samples cells per section and checks their cell/nucleus boundary polygons
directly (not just the pre-computed area columns): polygon validity,
nucleus-within-cell containment, and nucleus:cell area ratio. Flags
candidates for review; does not exclude anything (`04_quality_control/06_define_qc_thresholds_hierarchically.R`, `04_quality_control/07_apply_qc_filters_with_audit_trail.py`'s job).

Primary output: reports/qc/segmentation_review.html
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.segmentation_quality import build_segmentation_quality_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_path = project_root / "reports" / "qc" / "segmentation_review.html"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="03_assess_segmentation_quality",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_segmentation_quality_report(spatialdata_root, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   {summary['sections_processed']} section(s). "
        f"Median frac. nucleus-not-contained: {summary['median_fraction_nucleus_not_contained']:.4f}, "
        f"max: {summary['max_fraction_nucleus_not_contained']:.4f}. "
        f"{summary['sections_with_invalid_polygons']} section(s) with invalid polygons. "
        f"Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
