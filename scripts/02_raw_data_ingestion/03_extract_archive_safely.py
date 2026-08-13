#!/usr/bin/env python3
"""
`02_raw_data_ingestion/03_extract_archive_safely.py`

Extracts GSE300147_RAW.tar into data/staged/<GSM_accession>/, validating
every member path against traversal (absolute paths, ".." components)
before writing, and preserving original mtimes. Requires `02_raw_data_ingestion/02_verify_archive_checksums.py`'s
integrity check to have already passed (reads reports/integrity.tsv's
existence as a proxy -- does not re-verify checksums itself).

Primary output: data/staged/<sample_id>/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.ingest.extract import safe_extract


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    integrity_report = project_root / "reports" / "integrity.tsv"
    if not integrity_report.is_file():
        print(
            f"[ERROR] '{integrity_report}' not found. Run `02_raw_data_ingestion/02_verify_archive_checksums.py` (verify_archive_checksums) first.",
            file=sys.stderr,
        )
        return 1

    archive_path = project_root / "data" / "raw" / "GSE300147_RAW.tar"
    dest_root = project_root / "data" / "staged"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "02_raw_data_ingestion",
        script_name="03_extract_archive_safely",
        project_root=project_root,
        phase="02_raw_data_ingestion",
    )

    try:
        summary = safe_extract(archive_path, dest_root, project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(dest_root))
    logger.write(status="ok")
    print(
        f"[OK]   Extracted {summary['files_extracted']} file(s) across {summary['samples']} sample(s) to {dest_root}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
