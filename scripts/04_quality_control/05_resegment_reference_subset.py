#!/usr/bin/env python3
"""
`04_quality_control/05_resegment_reference_subset.py`

Independently resegments a representative subset of sections with an
alternative transcript-reassignment method (morphology-free nearest-
nucleus-centroid reassignment, 10x's own documented "nucleus expansion"
preset radius); used solely to test whether headline spatial results are
robust to segmentation choice. See
src/xenium_tcr_ecology/qc/resegmentation.py's module docstring
representative-subset selection, and scope documentation for this
milestone (deferred alongside 04.04 until Cell Type Annotation typing existed).

Primary output: data/objects/resegmented_subset/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.resegmentation import build_resegmentation_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="05_resegment_reference_subset",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_resegmentation_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['sections_processed']} section(s) resegmented. "
        f"Median total-count correlation: {summary['median_total_count_correlation']:.4f}, "
        f"median pseudobulk correlation: {summary['median_pseudobulk_correlation']:.4f}, "
        f"median fraction concordant (same cell): {summary['median_fraction_concordant_same_cell']:.4f}. "
        f"Wrote {summary['output_dir']}, {summary['summary_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
