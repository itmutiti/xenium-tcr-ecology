#!/usr/bin/env python3
"""
`06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R` helper -- _05_prepare_tme_substate_inputs.py

NOT its own numbered blueprint phase step: computes the marker-score
evidence 05_resolve_myeloid_and_stromal_substates.R needs for 5
compartments (Myeloid, Dendritic_cell, Fibroblast, Endothelial,
Perivascular_SmoothMuscle), since that R script cannot read
analysis_ready.h5ad directly. Invoked by the R script via a subprocess
call, not run standalone as a pipeline phase.

Primary output: data/derived/tme_substate_inputs.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.tme_substates import build_tme_substate_inputs


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
        script_name="_05_prepare_tme_substate_inputs",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_tme_substate_inputs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells, {summary['n_compartments']} compartment(s). "
        f"Cells per compartment: {summary['n_cells_per_compartment']}. "
        f"Wrote data/derived/tme_substate_inputs.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
