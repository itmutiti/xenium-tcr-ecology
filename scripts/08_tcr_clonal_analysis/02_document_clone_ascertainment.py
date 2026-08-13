#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`

Records how each probed clonotype was selected relative to each patient's
full repertoire and publishes this as an explicit boundary condition on
every downstream generalisability claim. The ascertainment criteria
are taken directly from the source paper's stated methods (abundance-
biased selection plus deliberate bystander/viral-reactive controls); this
project has no VDJ/TCR-contig data to independently reconstruct per-clone
repertoire rank, which is itself documented as the limitation. See
src/xenium_tcr_ecology/tcr/clone_ascertainment.py's module docstring.

Primary output: metadata/clone_ascertainment.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.clone_ascertainment import build_clone_ascertainment_record


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
        script_name="02_document_clone_ascertainment",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = build_clone_ascertainment_record(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_probes']} probe(s) documented. Source paper names "
        f"{summary['n_probes_source_paper_named_validated_clonotypes']} validated clonotypes across "
        f"{summary['n_patients_source_paper_named']} patients; this project covers "
        f"{summary['n_patients_this_project']} patients. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
