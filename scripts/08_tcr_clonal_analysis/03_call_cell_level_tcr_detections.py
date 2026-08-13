#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`

Calls clone detections from transcript counts with an explicit threshold
(any nonzero count, consistent with `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s own detection-rate
methodology), restricted to each probe's intended patient's own T cells
(`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`), and flags (does not resolve -- `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`'s job) multi-probe
ambiguous cells. See src/xenium_tcr_ecology/tcr/cell_calls.py's module
docstring.

Primary output: data/derived/tcr_cell_calls.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.cell_calls import build_tcr_cell_calls


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
        script_name="03_call_cell_level_tcr_detections",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_tcr_cell_calls(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_tcells']} T cell(s) evaluated against {summary['n_probes_used']} probe(s). "
        f"{summary['n_tcells_with_any_detection']} ({summary['fraction_tcells_with_detection']*100:.2f}%) "
        f"with a detection, {summary['n_tcells_multi_probe_ambiguous']} multi-probe ambiguous "
        f"({summary['n_tcells_likely_single_clone_tra_trb_pair']} of those look like a normal TRA+TRB pair, "
        f"{summary['n_tcells_multi_probe_ambiguous_excluding_likely_pairs']} ambiguous). "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
