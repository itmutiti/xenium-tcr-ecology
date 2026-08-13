#!/usr/bin/env python3
"""
`12_external_checkpoint_validation/03_decide_freeze_or_revise.py`

Applies a two-condition decision rule (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` module transfer
and `12_external_checkpoint_validation/02_quantify_directional_consistency.py` positive directional correlation) to decide freeze or
revise -- see
src/xenium_tcr_ecology/external_checkpoint/freeze_decision.py's module docstring
(including a note on why this rule was not declared strictly ahead of
the checks it gates).

Primary output: governance/freeze_decision.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.external_checkpoint.freeze_decision import build_freeze_decision


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "12_external_checkpoint_validation",
        script_name="03_decide_freeze_or_revise",
        project_root=project_root,
        phase="12_external_checkpoint_validation",
    )

    try:
        summary = build_freeze_decision(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Decision: {summary['decision'].upper()} (resulting version: {summary['resulting_version']}). "
        f"All programs transfer: {summary['all_programs_transfer']}. Overall Spearman rho: {summary['overall_spearman_rho']:.3f}. "
        f"Flagged states (caveat, not blocking): {summary['flagged_states']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
