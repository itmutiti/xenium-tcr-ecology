#!/usr/bin/env python3
"""
`02_raw_data_ingestion/05_standardise_sample_directory_layout.py`

Creates data/standardised/<section_id>/ for every staged sample, containing
symlinks (never copies) back to the files in data/staged/<GSM>/,
renamed to canonical role-based filenames -- so SpatialData Import onward reads a
predictable layout regardless of GEO's original, inconsistent embedded
timestamp filenames. section_id comes from metadata/sample_manifest.tsv
(`01_project_setup_and_governance/01_build_sample_manifest.py`), not re-derived here.

Primary output: data/standardised/<section_id>/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.ingest.standardize import standardize_layout


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    staged_root = project_root / "data" / "staged"
    standardised_root = project_root / "data" / "standardised"
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "02_raw_data_ingestion",
        script_name="05_standardise_sample_directory_layout",
        project_root=project_root,
        phase="02_raw_data_ingestion",
    )

    try:
        summary = standardize_layout(
            staged_root, standardised_root, sample_manifest_path, project_root
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(standardised_root))
    logger.write(status="ok")
    print(
        f"[OK]   Standardised {summary['sections_standardised']} section(s) into {standardised_root}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
