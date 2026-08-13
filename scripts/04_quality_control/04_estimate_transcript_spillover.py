#!/usr/bin/env python3
"""
`04_quality_control/04_estimate_transcript_spillover.py`

Flags cells near segmentation boundaries adjacent to a different predicted
cell type; estimates a per-cell spillover-risk score using cell-
boundary-polygon distance-to-neighbour and neighbour-identity weighting
against `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s final lineage calls.

This milestone was deferred until after Cell Type Annotation: the blueprint's own spec
requires a predicted cell type per cell, which did not exist earlier in the
pipeline. Not a registered hypothesis test (governance/analysis_registry.tsv
covers hypothesis-bearing/formal validation analyses only; this is a QC/
engineering step, same category as 04.00-04.03).

Primary output: data/derived/spillover_risk.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.qc.transcript_spillover import build_spillover_risk_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "04_quality_control",
        script_name="04_estimate_transcript_spillover",
        project_root=project_root,
        phase="04_quality_control",
    )

    try:
        summary = build_spillover_risk_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} cell(s) across {summary['n_sections']} section(s). "
        f"Mean spillover risk score: {summary['mean_spillover_risk_score']:.4f}. "
        f"{summary['fraction_boundary_adjacent_to_different_type']*100:.2f}% boundary-adjacent "
        f"to a different predicted type. Wrote data/derived/spillover_risk.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
