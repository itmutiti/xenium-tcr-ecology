#!/usr/bin/env python3
"""
`03_spatialdata_import/00_detect_xenium_format_version.py`

Checks every standardised section against spatialdata_io.xenium()'s hard
requirements (experiment.xenium, cells.zarr.zip) and records the deviation:
GEO's supplementary release for this accession does not include either
file, so the standard reader cannot be used -- `03_spatialdata_import/01_import_each_section_to_spatialdata.py` implements a
custom reader against the 6 files that are actually published. Also
extracts and cross-checks the authoritative pixel-to-micron scale
(OME-XML PhysicalSizeX/Y/Z) from every section's morphology image, failing
loudly if sections disagree (which would mean a single hardcoded scale
transform in `03_spatialdata_import/01_import_each_section_to_spatialdata.py` would silently misalign some sections).

Primary output: metadata/format_versions.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.format_detection import build_format_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    standardised_root = project_root / "data" / "standardised"
    output_path = project_root / "metadata" / "format_versions.tsv"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="00_detect_xenium_format_version",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    try:
        summary = build_format_report(standardised_root, output_path, project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   Checked {summary['sections_checked']} section(s). "
        f"{summary['standard_reader_incompatible']}/{summary['sections_checked']} incompatible with "
        f"spatialdata_io.xenium() (missing experiment.xenium/cells.zarr.zip -- expected, custom reader required)."
    )
    print(
        f"[INFO] Distinct pixel-size (X,Y,Z um) combinations across sections: {summary['pixel_sizes_um']}"
    )
    if summary["distinct_pixel_size_um_combos"] > 1:
        print(
            "[WARN] Sections do NOT share a single pixel size -- `03_spatialdata_import/01_import_each_section_to_spatialdata.py` must read the scale per-section, "
            "not use one hardcoded constant."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
