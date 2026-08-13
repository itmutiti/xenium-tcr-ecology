#!/usr/bin/env python3
"""
`06_cell_type_annotation/04_resolve_t_cell_substates.R` helper -- _04_prepare_t_cell_substate_inputs.py

NOT its own numbered blueprint phase step: computes the marker evidence
(Treg score, CD4/CD8A expression, reused `05_preprocessing_and_normalisation/03_calculate_program_scores.py` program scores, Phase
6.02 lineage argmax, `06_cell_type_annotation/03_map_external_scrna_reference.py` reference-transfer prediction) that
04_resolve_t_cell_substates.R needs, since that R script cannot read
analysis_ready.h5ad directly. Invoked by the R script via a subprocess
call, not run standalone as a pipeline phase.

Primary output: data/derived/t_cell_substate_inputs.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.t_cell_substates import build_t_cell_substate_inputs


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
        script_name="_04_prepare_t_cell_substate_inputs",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_t_cell_substate_inputs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells, {summary['n_t_or_nk_lineage_cells']:,} T/NK-lineage cells "
        f"(`06_cell_type_annotation/02_score_major_lineages.py` argmax). Wrote data/derived/t_cell_substate_inputs.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
