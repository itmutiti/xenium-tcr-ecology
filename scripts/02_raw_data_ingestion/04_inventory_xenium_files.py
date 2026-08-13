#!/usr/bin/env python3
"""
`02_raw_data_ingestion/04_inventory_xenium_files.py`

Scans data/staged/<GSM>/ for the 6 mandatory Xenium output files per sample
(cell_boundaries, cell_feature_matrix, cells, morphology, nucleus_boundaries,
transcripts -- roles verified against this accession's filelist.txt,
not assumed) plus the optional 7th raw-microscopy TIFF, and fails loudly if
any sample is missing a mandatory file rather than silently proceeding.

Primary output: metadata/xenium_file_inventory.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.ingest.xenium_inventory import build_xenium_file_inventory


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    staged_root = project_root / "data" / "staged"
    output_path = project_root / "metadata" / "xenium_file_inventory.tsv"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "02_raw_data_ingestion",
        script_name="04_inventory_xenium_files",
        project_root=project_root,
        phase="02_raw_data_ingestion",
    )

    try:
        summary = build_xenium_file_inventory(staged_root, output_path, project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   {summary['samples_complete']}/{summary['samples_scanned']} sample(s) have all mandatory files. Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
