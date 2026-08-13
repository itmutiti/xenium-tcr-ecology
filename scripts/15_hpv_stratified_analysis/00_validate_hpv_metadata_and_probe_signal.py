#!/usr/bin/env python3
"""
`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`

Cross-checks clinical p16 IHC labels against Xenium HPV16 E6/E7
oncogene probe signal, flagging discordant or clinically-unverifiable
patients BEFORE `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py` fixes the primary HPV contrast. See
src/xenium_tcr_ecology/hpv/hpv_validation.py's module docstring.

Primary output: metadata/hpv_status_validated.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.hpv.hpv_validation import build_hpv_status_validation


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
        script_name="00_validate_hpv_metadata_and_probe_signal",
        project_root=project_root,
        phase="15_hpv_stratified_analysis",
    )

    try:
        summary = build_hpv_status_validation(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_patients']} patient(s) validated: {summary['n_confirmed_positive']} confirmed positive, "
        f"{summary['n_confirmed_negative_or_no_verification']} confirmed/presumed negative, "
        f"{summary['n_discordant']} DISCORDANT, {summary['n_probe_positive_clinically_untested']} probe-positive-but-untested, "
        f"{summary['n_presumed_negative_unverifiable']} presumed-negative-unverifiable. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
