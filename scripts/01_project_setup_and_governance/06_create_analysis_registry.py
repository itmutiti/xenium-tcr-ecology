#!/usr/bin/env python3
"""
`01_project_setup_and_governance/06_create_analysis_registry.py`

Compiles config/governance/analysis_registry_input.yaml into
governance/analysis_registry.tsv: one row per planned hypothesis-bearing or
formal validation analysis (not every phase script -- most are QC/
engineering steps with no hypothesis to register), each with its unit of
analysis, exclusion criteria, primary endpoint, and multiplicity family
pre-committed before any of these analyses has actually been run.

Enforces the HPV single-contrast cap (config/reproducibility_policy.yaml
gates.hpv_single_contrast) at registration time, not just at HPV-Stratified Analysis.

Primary output: governance/analysis_registry.tsv
"""

from __future__ import annotations

import sys
from datetime import date

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.governance.registry import compile_analysis_registry
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root


def main() -> int:
    parser = base_parser(__doc__)
    parser.add_argument("--registered-by", default="Irvine Tatenda Mutiti")
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    input_path = project_root / "config" / "governance" / "analysis_registry_input.yaml"
    output_path = project_root / "governance" / "analysis_registry.tsv"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "01_project_setup_and_governance",
        script_name="06_create_analysis_registry",
        project_root=project_root,
        phase="01_project_setup_and_governance",
    )

    try:
        summary = compile_analysis_registry(
            input_path,
            output_path,
            project_root,
            registered_by=args.registered_by,
            registered_date=date.today().isoformat(),
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(f"[OK]   Registered {summary['analyses_registered']} analys(es) to {output_path}")
    print(
        f"[INFO] {summary['hpv_primary_contrasts_reserved']} HPV primary contrast slot(s) reserved "
        f"(cap: 2); {summary['distinct_primary_families']} distinct primary multiplicity famil(y/ies)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
