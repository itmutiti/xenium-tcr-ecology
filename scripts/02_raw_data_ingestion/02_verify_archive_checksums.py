#!/usr/bin/env python3
"""
`02_raw_data_ingestion/02_verify_archive_checksums.py`

Computes this project's own SHA-256 of the downloaded archive (GEO
publishes no per-file checksum for this accession, so there is nothing
external to verify against) and confirms the archive's internal file
listing (read without extracting) exactly matches the filelist.txt
manifest recorded in `02_raw_data_ingestion/00_query_geo_accession.py`'s geo_snapshot.json -- before `02_raw_data_ingestion/03_extract_archive_safely.py`
spends time and disk extracting a possibly-truncated download.

Primary output: data/raw/SHA256SUMS; reports/integrity.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.ingest.integrity import verify_archive_and_checksum


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    archive_path = project_root / "data" / "raw" / "GSE300147_RAW.tar"
    geo_snapshot_path = project_root / "metadata" / "geo_snapshot.json"
    sha256sums_path = project_root / "data" / "raw" / "SHA256SUMS"
    integrity_report_path = project_root / "reports" / "integrity.tsv"

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "02_raw_data_ingestion",
        script_name="02_verify_archive_checksums",
        project_root=project_root,
        phase="02_raw_data_ingestion",
    )

    try:
        summary = verify_archive_and_checksum(
            archive_path, geo_snapshot_path, sha256sums_path, integrity_report_path, project_root
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(f"[OK]   Archive SHA-256: {summary['archive_sha256']}")
    print(
        f"[OK]   {summary['files_verified']}/{summary['files_expected']} files verified against filelist.txt. "
        f"Wrote {sha256sums_path} and {integrity_report_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
