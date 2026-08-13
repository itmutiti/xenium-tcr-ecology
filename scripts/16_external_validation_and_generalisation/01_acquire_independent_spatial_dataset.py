#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py`

Downloads (if not already present) and verifies the independent Xenium
spatial dataset (Janesick et al. 2023, Nat Commun -- see `data/external/
spatial/Xenium_Janesick_BreastCancer_Rep1/README.md` for the full
acquisition provenance, citation and license). Public, unauthenticated
10x Genomics CDN download, no manual step required. See
src/xenium_tcr_ecology/validation/spatial_dataset_acquisition.py's
module docstring.

Primary output: data/external/spatial/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.spatial_dataset_acquisition import (
    build_spatial_dataset_acquisition_summary,
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
        script_name="01_acquire_independent_spatial_dataset",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_spatial_dataset_acquisition_summary(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Independent spatial dataset verified: {summary['n_cells']} cells, {summary['n_clusters']} clusters, "
        f"{summary['n_files_verified']} file(s) checksum-verified. {summary['dataset_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
