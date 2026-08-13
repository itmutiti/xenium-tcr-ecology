#!/usr/bin/env python3
"""
`05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`

Consolidates analysis_ready.h5ad (`05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`) with program_scores.parquet
(`05_preprocessing_and_normalisation/03_calculate_program_scores.py`) into one object, freezes it plus supporting diagnostic tables
into a versioned release directory, and records a SHA256 manifest of
exactly what was used.

Makes no new scientific decisions -- `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` already selected the
primary normalisation layer, per its documented evidence-based comparison; this
is a mechanical, verifiable freeze.

Primary output: data/releases/v1_primary_analysis/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.preprocess.release_freeze import (
    RELEASE_NAME,
    freeze_primary_analysis_matrix,
)


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    release_dir = project_root / "data" / "releases" / RELEASE_NAME
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "05_preprocessing_and_normalisation",
        script_name="05_create_primary_analysis_matrix",
        project_root=project_root,
        phase="05_preprocessing_and_normalisation",
    )

    try:
        summary = freeze_primary_analysis_matrix(project_root, release_dir)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells x {summary['n_genes']} genes, "
        f"primary layer='{summary['primary_normalization_layer']}'. "
        f"{summary['n_files']} file(s) frozen with SHA256 checksums. "
        f"Matrix hash: {summary['matrix_hash'][:16]}... "
        f"Wrote {summary['release_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
