#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`

Classifies every `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`-evaluated T cell into unassigned, singlet,
probable_multiplet, or low_confidence, combining `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s detection
calls with `08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`'s per-probe empirical false-positive-rate
estimates. Cross-checks the result against `04_quality_control/04_estimate_transcript_spillover.py`'s independent
transcript-spillover-risk score. See
src/xenium_tcr_ecology/tcr/resolve_calls.py's module docstring.

Primary output: data/derived/tcr_resolved_calls.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.resolve_calls import build_resolved_calls


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
        script_name="06_resolve_multiclonal_and_ambiguous_cells",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_resolved_calls(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} T cell(s) resolved: {summary['resolution_counts']}. "
        f"Mean spillover_risk_score by resolution: {summary['mean_spillover_risk_by_resolution']}. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
