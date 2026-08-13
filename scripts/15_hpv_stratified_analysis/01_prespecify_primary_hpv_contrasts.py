#!/usr/bin/env python3
"""
`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`

Fixes the single primary, molecularly-validated HPV-positive vs
HPV-negative contrast for the entire thesis, using `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`'s
`validated_hpv_status` (not the raw clinical label), before any Phase
15 hypothesis test runs. See src/xenium_tcr_ecology/hpv/primary_
contrast.py's module docstring.

Primary output: governance/hpv_primary_contrasts.yaml
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.hpv.primary_contrast import build_primary_hpv_contrasts


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
        script_name="01_prespecify_primary_hpv_contrasts",
        project_root=project_root,
        phase="15_hpv_stratified_analysis",
    )

    try:
        summary = build_primary_hpv_contrasts(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Primary HPV contrast '{summary['contrast_name']}' fixed: "
        f"n={summary['n_positive']} positive vs n={summary['n_negative']} negative "
        f"({summary['n_excluded']} excluded, discordant/unverifiable). Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
