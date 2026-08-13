#!/usr/bin/env python3
"""
`14_spatial_interactions_and_barriers/05_prioritise_testable_interactions.py`

Ranks `14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`'s candidate spatial interactions by effect
size, spatial specificity, cross-patient consistency, panel completeness
and a qualitative external-support tier. Ranks spatial associations for
prioritised follow-up testing -- makes no claim about which pathway
underlies any association. See
src/xenium_tcr_ecology/interactions/prioritisation.py's module docstring
.

Primary output: results/interaction_priority_table.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.interactions.prioritisation import build_interaction_priority_table


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "14_spatial_interactions_and_barriers",
        script_name="05_prioritise_testable_interactions",
        project_root=project_root,
        phase="14_spatial_interactions_and_barriers",
    )

    try:
        summary = build_interaction_priority_table(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_candidate_interactions']} candidate interaction(s) ranked. Top priority: "
        f"{summary['top_priority_sender_receiver_pair_id']} / {summary['top_priority_lr_pair_id']}. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
