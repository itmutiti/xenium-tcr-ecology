#!/usr/bin/env python3
"""
`03_spatialdata_import/02_validate_coordinate_systems.py`

Extends `03_spatialdata_import/01_import_each_section_to_spatialdata.py`'s single-section manual alignment check (DAPI intensity
at cell centroids vs. random background) into an automated,
per-section check across every imported SpatialData store: confirms cell
counts match boundary counts, every element sits under one unified "global"
coordinate system, and image/points/shapes are physically aligned -- not
just structurally present. Fails loudly (non-zero exit) if any section
fails, rather than silently producing a report nobody reads.

Primary output: reports/coordinate_validation/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.coordinate_validation import build_coordinate_validation_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    format_versions_path = project_root / "metadata" / "format_versions.tsv"
    output_dir = project_root / "reports" / "coordinate_validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "validation_report.tsv"

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="02_validate_coordinate_systems",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    try:
        summary = build_coordinate_validation_report(
            spatialdata_root, format_versions_path, output_path, project_root
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   {summary['sections_passed']}/{summary['sections_validated']} section(s) passed coordinate validation. "
        f"Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
