#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`

Downloads (if not already present) and verifies the HNSCC reference
datasets: GSE139324 (Cillo et al. 2020, Immunity -- a second independent
HNSCC scRNA-seq reference, distinct from GSE103322 already used in
External Checkpoint Validation) and TCGA-HNSC (bulk RNA-seq +
clinical/survival, via UCSC Xena). Both public, unauthenticated
downloads, no manual step required. See
`data/external/scrna/GSE139324/README.md` and
`data/external/bulk/TCGA-HNSC/README.md` for the full acquisition
provenance.

Primary output: data/external/scrna/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.scrna_reference_acquisition import (
    build_scrna_reference_acquisition_summary,
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
        logs_dir=project_root / "results" / "logs" / "16_external_validation_and_generalisation",
        script_name="02_acquire_hnscc_scrna_references",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_scrna_reference_acquisition_summary(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   HNSCC reference datasets verified: GSE139324 "
        f"{summary['gse139324_n_til_samples']} TIL samples ({summary['gse139324_n_files_verified']} file(s) checksummed), "
        f"TCGA-HNSC {summary['tcga_n_samples']} samples ({summary['tcga_n_files_verified']} file(s) checksummed)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
