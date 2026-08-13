#!/usr/bin/env python3
"""
`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` helper -- _02_compute_normalization_benchmark_metrics.py

NOT its own numbered blueprint phase step: computes the replicate-stability
and negative-control-probe-behaviour metrics that
02_evaluate_normalisation_strategies.R needs, since that R script cannot
read analysis_ready.h5ad directly (no R HDF5/AnnData reader available). Invoked by the R script
via a subprocess call, not run standalone as a pipeline phase.

Primary output: reports/preprocess/normalisation_benchmark_replicate_stability.parquet,
reports/preprocess/normalisation_benchmark_technical_noise.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.preprocess.normalization_benchmark import (
    build_normalization_benchmark_summary,
)


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    output_path = project_root / "reports" / "preprocess" / "normalisation_benchmark"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "05_preprocessing_and_normalisation",
        script_name="_02_compute_normalization_benchmark_metrics",
        project_root=project_root,
        phase="05_preprocessing_and_normalisation",
    )

    try:
        summary = build_normalization_benchmark_summary(project_root, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_replicate_pairs']} pair(s). "
        f"Median replicate r by method: {summary['median_replicate_r_by_method']}. "
        f"Abs. technical-noise rho by method: {summary['abs_technical_noise_rho_by_method']}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
