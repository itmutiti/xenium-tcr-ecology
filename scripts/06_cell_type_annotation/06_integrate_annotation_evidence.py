#!/usr/bin/env python3
"""
`06_cell_type_annotation/06_integrate_annotation_evidence.py`

Combines marker (`06_cell_type_annotation/02_score_major_lineages.py`), cluster (`06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`), and spatial-neighbourhood
evidence into a final major-lineage call with a composite [0,1] confidence
score and an explicit ambiguity flag; folds in `06_cell_type_annotation/04_resolve_t_cell_substates.R`, `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R`'s substate
calls, cleared for ambiguous cells per the blueprint's "fine subtypes are
withheld when the panel cannot support them." See
src/xenium_tcr_ecology/annotation/integrate_evidence.py's module docstring
for the full methodology.

Primary output: data/derived/final_cell_annotations.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.integrate_evidence import build_final_annotations_report


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
        script_name="06_integrate_annotation_evidence",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_final_annotations_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells. {summary['n_ambiguous']:,} ambiguous "
        f"({100 * summary['fraction_ambiguous']:.2f}%). Mean confidence: {summary['mean_confidence']:.3f}. "
        f"{summary['n_with_substate']:,} cells retain a substate call. "
        f"Final lineage counts: {summary['final_lineage_counts']}. "
        f"Wrote data/derived/final_cell_annotations.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
