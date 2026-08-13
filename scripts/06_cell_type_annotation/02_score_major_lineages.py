#!/usr/bin/env python3
"""
`06_cell_type_annotation/02_score_major_lineages.py`

Assigns per-cell scores for each of `06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s major_lineage identities
using scanpy's score_genes on `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`'s primary normalisation layer.
Scores only, not final cell-type calls -- integrating this with clustering,
reference-mapping and spatial evidence is `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s job.

Primary output: data/derived/lineage_scores.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.lineage_scores import build_lineage_scores_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "06_cell_type_annotation",
        script_name="02_score_major_lineages",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_lineage_scores_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells, {summary['n_lineages']} lineage(s). "
        f"Argmax lineage counts: {summary['argmax_lineage_counts']}. "
        f"Wrote data/derived/lineage_scores.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
