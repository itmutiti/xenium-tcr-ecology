#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`

Ensures patient-specific probes are only evaluated in their intended
specimens and detects leakage or naming conflicts: for each probe,
determines its most likely intended patient from T-cell detection
evidence (Fisher's exact test, FDR-corrected across all probes), since
`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py` found every probe is physically present on 3-4 patients'
panels by manufacturing-batch structure, not by design intent. See
src/xenium_tcr_ecology/tcr/patient_mapping.py's module docstring.

Primary output: reports/tcr/patient_probe_audit.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.patient_mapping import build_patient_probe_audit


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
        script_name="01_map_tcr_probes_to_patients",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_patient_probe_audit(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_probes_audited']} probe(s) audited using {summary['n_tcells_used']} T cell(s). "
        f"{summary['n_probes_with_identified_patient']} with a statistically identified intended patient "
        f"(FDR alpha={summary['fdr_alpha']}), {summary['n_probes_no_significant_specificity']} with no "
        f"significant patient specificity. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
