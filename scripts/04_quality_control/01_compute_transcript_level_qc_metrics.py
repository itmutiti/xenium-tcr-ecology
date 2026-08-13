#!/usr/bin/env python3
"""
`04_quality_control/01_compute_transcript_level_qc_metrics.py`

Profiles Q-values, unassigned-transcript rate, nuclear overlap, and
negative-control/unassigned-codeword burden per section, reading from each
section's `transcripts` Points element (not `table` -- the cell x
gene matrix excludes these feature types but the transcripts element does
not).

Primary output: data/derived/transcript_qc_metrics.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.transcript_metrics import build_transcript_qc_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_path = project_root / "data" / "derived" / "transcript_qc_metrics.parquet"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="01_compute_transcript_level_qc_metrics",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_transcript_qc_report(spatialdata_root, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   {summary['sections_processed']} section(s), {summary['total_transcripts']:,} transcripts total. "
        f"Median fraction QV<20: {summary['median_fraction_qv_below_20']:.3f}, "
        f"median fraction overlapping nucleus: {summary['median_fraction_overlaps_nucleus']:.3f}. Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
