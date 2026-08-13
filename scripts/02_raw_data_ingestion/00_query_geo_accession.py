#!/usr/bin/env python3
"""
`02_raw_data_ingestion/00_query_geo_accession.py`

Fetches GSE300147's filelist.txt from the GEO FTP suppl directory
(public, unauthenticated) and records a timestamped snapshot: archive size,
per-file sizes/types, retrieval date, and the GSM sample list re-derived
from file naming -- independent of, and cross-checkable against, Phase
1.01's manifest (built from individual GSM SOFT records, a different GEO
endpoint).

Primary output: metadata/geo_snapshot.json
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.ingest.geo_query import build_geo_snapshot
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root

ACCESSION = "GSE300147"


def main() -> int:
    parser = base_parser(__doc__)
    parser.add_argument("--accession", default=ACCESSION)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    output_path = project_root / "metadata" / "geo_snapshot.json"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "02_raw_data_ingestion",
        script_name="00_query_geo_accession",
        project_root=project_root,
        phase="02_raw_data_ingestion",
    )

    try:
        summary = build_geo_snapshot(args.accession, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    gb = summary["archive_size_bytes"] / 1e9
    print(
        f"[OK]   {args.accession}: archive {gb:.2f} GB, {summary['sample_count']} samples, "
        f"{summary['file_count']} files. Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
