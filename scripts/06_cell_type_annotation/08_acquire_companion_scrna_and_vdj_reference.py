#!/usr/bin/env python3
"""
`06_cell_type_annotation/08_acquire_companion_scrna_and_vdj_reference.py`

Downloads (if not already present) and verifies GSE287301 -- McCord et
al. 2026's own companion scRNA-seq dataset (`data/external/GSE287301/
README.md`): the aggregated gene-expression matrix and the 16 per-sample
paired scTCR-seq VDJ archives (`data/external/GSE287301/vdj/README.md`).
Public, unauthenticated NCBI GEO download, no manual step required.

Numbered as a Cell Type Annotation milestone (its first consumer,
`06_cell_type_annotation/03_map_external_scrna_reference.py`) rather than
renumbered into an earlier position, to avoid disturbing any existing
rule's identity -- `06_cell_type_annotation/03_...` and
`08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_vdj_
ground_truth.py` both depend on this milestone's sentinel directly (see
their own Snakemake rules' `input:`), not on numeric proximity.

Primary output: data/external/GSE287301/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.companion_reference_acquisition import (
    build_companion_reference_acquisition_summary,
)


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "06_cell_type_annotation",
        script_name="08_acquire_companion_scrna_and_vdj_reference",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_companion_reference_acquisition_summary(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   GSE287301 companion reference verified: gene-expression matrix at "
        f"{summary['gex_dir']}, {summary['vdj_pools_acquired']}/{summary['vdj_pools_expected']} "
        f"VDJ pools at {summary['vdj_dir']}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
