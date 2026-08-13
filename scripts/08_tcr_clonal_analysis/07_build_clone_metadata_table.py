#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/07_build_clone_metadata_table.py`

Summarises clone size, patient, section support, phenotype composition
and replicate recurrence for every distinct clone (`08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`'s `singlet`/
`low_confidence` cells, grouped by detected-probe identity). See
src/xenium_tcr_ecology/tcr/clone_metadata.py's module docstring
internal-consistency check (a clone must not span more than one
patient) and the replicate-recurrence check reusing `04_quality_control/08_assess_replicate_concordance.R`'s
technical-replicate section pairs.

Primary output: data/derived/clone_metadata.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.clone_metadata import build_clone_metadata


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "08_tcr_clonal_analysis",
        script_name="07_build_clone_metadata_table",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_clone_metadata(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clones']} clone(s), {summary['n_cells_in_clones']} cell(s). "
        f"Median clone size {summary['median_clone_size']}, max {summary['max_clone_size']}. "
        f"{summary['n_clones_detected_in_both_replicates']}/{summary['n_clones_in_replicate_patients']} "
        f"replicate-patient clones detected in both runs. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
