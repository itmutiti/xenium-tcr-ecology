#!/usr/bin/env python3
"""
`15_hpv_stratified_analysis/06_prepare_hpv_claim_strength_table.py`

Grades every HPV-related conclusion from `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`-15.05 as
`supported`, `exploratory`, or `unsuitable_for_inference` -- the final
synthesis this deliberately small (n=4 vs n=4) design requires.
See src/xenium_tcr_ecology/hpv/claim_strength.py's module docstring.

Primary output: results/hpv_claim_strength.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.hpv.claim_strength import build_hpv_claim_strength_table


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "15_hpv_stratified_analysis",
        script_name="06_prepare_hpv_claim_strength_table",
        project_root=project_root,
        phase="15_hpv_stratified_analysis",
    )

    try:
        summary = build_hpv_claim_strength_table(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_claims']} HPV claim(s) graded: {summary['n_supported']} supported, "
        f"{summary['n_exploratory']} exploratory, {summary['n_unsuitable_for_inference']} unsuitable for inference. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
