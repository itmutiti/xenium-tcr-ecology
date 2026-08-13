#!/usr/bin/env python3
"""
`04_quality_control/02_detect_spatial_qc_artifacts.py`

Flags candidate local-decoding-failure / striping artefacts at the
field-of-view (FOV) level -- confirmed against data as the natural
spatial tiling unit for this platform (27 FOVs/section, grid-named e.g.
'M6'). Uses a robust (median-absolute-deviation based) outlier score per
FOV within each section, rather than excluding anything: exclusion with a
documented rationale is `04_quality_control/06_define_qc_thresholds_hierarchically.R`, `04_quality_control/07_apply_qc_filters_with_audit_trail.py`'s job.

Primary output: reports/qc/spatial_artifact_masks/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.spatial_artifacts import build_spatial_artifact_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_dir = project_root / "reports" / "qc" / "spatial_artifact_masks"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="02_detect_spatial_qc_artifacts",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_spatial_artifact_report(spatialdata_root, output_dir, project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_dir))
    logger.write(status="ok")
    print(
        f"[OK]   {summary['total_fovs']} FOV(s) across {summary['sections_processed']} section(s). "
        f"{summary['total_flagged_fovs']} flagged as artefact candidates, across "
        f"{summary['sections_with_flagged_fovs']} section(s). Wrote {output_dir}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
