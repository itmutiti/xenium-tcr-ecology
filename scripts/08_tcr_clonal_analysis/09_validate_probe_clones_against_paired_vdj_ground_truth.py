#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_vdj_ground_truth.py`

Compares Xenium CDR3-probe patient assignments
(`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`) against paired
scTCR-seq VDJ data from the same patients (GSE287301): whether a probed
CDR3 sequence is detected in that patient's own VDJ repertoire, and
whether Xenium's spatial detection rate tracks clonal abundance rank.
Added after acquiring VDJ data not available at the time of the original
TCR pipeline design; see `docs/analysis_amendments.md`.

Method: `src/xenium_tcr_ecology/tcr/vdj_ground_truth_validation.py`.

Primary output: data/derived/probe_vdj_ground_truth_comparison.parquet,
reports/tcr/probe_vdj_ground_truth_validation.pdf
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.vdj_ground_truth_validation import (
    compare_probe_detections_to_vdj_ground_truth,
)


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
        script_name="09_validate_probe_clones_against_paired_vdj_ground_truth",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = compare_probe_detections_to_vdj_ground_truth(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_found_in_real_vdj_ground_truth']}/{summary['n_probes_with_identified_patient']} "
        f"({summary['fraction_found']:.1%}) patient-identified probes confirmed in the paired VDJ data."
    )
    print(
        f"[OK]   Paired VDJ data: {summary['n_real_vdj_ground_truth_cells_total']} cells, {summary['n_real_vdj_ground_truth_patients']} patients."
    )
    print(
        f"[OK]   Xenium detection rate vs. VDJ clonal-abundance rank, Spearman correlation: {summary['xenium_detection_vs_vdj_rank_spearman']}"
    )
    print(
        f"[OK]   Wrote {summary['output_path']}, {summary['ground_truth_path']}, {summary['report_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
