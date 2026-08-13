#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/05_screen_cdr3_cross_patient_similarity.py`

Screens all same-chain pairs of probed CDR3 sequences (`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`'s
registry) for high sequence similarity (Levenshtein edit distance <= 2)
that could indicate probe cross-reactivity, flagging cross-patient hits
specifically -- an expected phenomenon for public/quasi-public
viral-reactive motifs (the source paper's Figure 3 deliberately
includes VDJdb-matched bystander clones), not automatically a design
error. See src/xenium_tcr_ecology/tcr/cdr3_similarity.py's module
docstring.

Primary output: reports/tcr/cdr3_similarity_screen.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.cdr3_similarity import build_cdr3_similarity_screen


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
        script_name="05_screen_cdr3_cross_patient_similarity",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_cdr3_similarity_screen(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_probes']} probe(s), {summary['n_same_chain_pairs_screened']} same-chain pair(s) "
        f"screened (edit distance <= {summary['max_similarity_edit_distance']}). "
        f"{summary['n_similar_pairs_flagged']} similar pair(s) flagged, "
        f"{summary['n_cross_patient_similar_pairs']} cross-patient "
        f"({summary['n_probes_involved_in_cross_patient_similarity']} distinct probes involved). "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
