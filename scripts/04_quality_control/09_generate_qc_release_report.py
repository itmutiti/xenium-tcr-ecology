#!/usr/bin/env python3
"""
`04_quality_control/09_generate_qc_release_report.py`

Aggregates `04_quality_control/00_compute_cell_level_qc_metrics.py`-4.03 and 4.06-4.08 into a section-by-section QC report
and a single go/no-go decision. The decision is a CONDITIONAL GO: Phase
4.04/4.05 are deferred,
so Preprocessing and Normalisation/6 are unblocked but Clone Spatial Descriptors+/13 tumour-engagement claims remain
blocked until a follow-up release amends this decision.

Primary output: reports/qc/QC_release_report.html
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.qc_release_report import build_qc_release_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    output_path = project_root / "reports" / "qc" / "QC_release_report.html"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="09_generate_qc_release_report",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_qc_release_report(project_root, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Decision: {summary['status']}. {summary['sections_processed']} section(s), "
        f"{summary['n_cells_raw_total']:,} -> {summary['n_cells_post_qc_total']:,} cells. "
        f"Outstanding: {len(summary['outstanding_workstreams'])} workstream(s). "
        f"Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
