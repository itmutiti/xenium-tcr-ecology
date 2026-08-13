#!/usr/bin/env python3
"""
`04_quality_control/07_apply_qc_filters_with_audit_trail.py`

Applies `04_quality_control/06_define_qc_thresholds_hierarchically.R`'s active QC threshold profile plus `04_quality_control/02_detect_spatial_qc_artifacts.py`'s
FOV-artifact flags to every cell, records a reason code for every excluded
cell, and writes a filtered analysis object. The untouched `03_spatialdata_import/05_build_combined_analysis_object.py`
combined object remains on disk as the immutable source of truth.

Transcript-level QV filtering was considered and explicitly NOT applied to
the primary count matrix (human-approved scope decision).

Primary output: data/objects/qc_filtered.h5ad; data/derived/exclusion_log.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.apply_filters import build_qc_filter_report


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
        script_name="07_apply_qc_filters_with_audit_trail",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_qc_filter_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Profile '{summary['active_profile']}': "
        f"{summary['n_cells_excluded']} / {summary['n_cells_total']} cells excluded "
        f"({100 * summary['fraction_excluded']:.2f}%), {summary['n_cells_retained']} retained. "
        f"Reason breakdown: {summary['reason_counts']}. "
        f"Wrote data/objects/qc_filtered.h5ad and data/derived/exclusion_log.tsv"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
