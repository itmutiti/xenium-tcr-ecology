#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`

Separates patient-specific CDR3 probes from conventional T-cell genes
using `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s feature metadata and the probe naming convention
(date/batch prefix, CDR3 amino acid sequence, TCR chain), then structures
each probe with the section(s)/patient(s) whose panel physically includes
it. See src/xenium_tcr_ecology/tcr/probe_registry.py's module docstring
for scope (this builds the registry; auditing it for leakage is Phase
8.01's job, not pre-empted here).

Primary output: metadata/tcr_probe_registry.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.probe_registry import build_tcr_probe_registry


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
        script_name="00_identify_tcr_cdr3_probe_features",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_tcr_probe_registry(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_probes']} probe(s): {summary['n_tra_probes']} TRA, {summary['n_trb_probes']} TRB, "
        f"{summary['n_distinct_date_batch_prefixes']} distinct date/batch prefixes. "
        f"{summary['n_probes_single_patient']} single-patient, {summary['n_probes_multi_patient']} multi-patient, "
        f"{summary['n_probes_zero_patients']} zero-patient probes. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
