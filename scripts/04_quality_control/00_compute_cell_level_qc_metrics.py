#!/usr/bin/env python3
"""
`04_quality_control/00_compute_cell_level_qc_metrics.py`

Computes per-cell QC metrics from the `03_spatialdata_import/05_build_combined_analysis_object.py` combined object:
transcript counts (cross-checked against the expression matrix, not just
trusted from cells.parquet), detected genes, cell/nucleus area, control
probe/codeword ratios, and transcript density.

Primary output: data/derived/cell_qc_metrics.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.cell_metrics import build_cell_qc_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    combined_h5ad_path = project_root / "data" / "objects" / "hnscc_xenium_combined.h5ad"
    output_path = project_root / "data" / "derived" / "cell_qc_metrics.parquet"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="00_compute_cell_level_qc_metrics",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_cell_qc_report(combined_h5ad_path, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} cells across {summary['n_sections']} sections. "
        f"Median transcripts/cell: {summary['median_transcript_counts']:.0f}, "
        f"median genes/cell: {summary['median_genes_detected']:.0f}. Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
