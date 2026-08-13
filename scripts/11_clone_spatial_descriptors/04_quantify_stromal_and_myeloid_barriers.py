#!/usr/bin/env python3
"""
`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`

Computes clone-specific fibroblast, vascular and suppressive-
myeloid barrier composition along the shortest graph path to the
nearest tumour cell, each with a calibrated (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`-style
constrained-permutation) empirical p-value -- see
src/xenium_tcr_ecology/clone_ecology/barrier_metrics.py's module
docstring.

Primary output: data/derived/clone_barrier_metrics.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.barrier_metrics import build_clone_barrier_metrics


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "11_clone_spatial_descriptors",
        script_name="04_quantify_stromal_and_myeloid_barriers",
        project_root=project_root,
        phase="11_clone_spatial_descriptors",
    )

    try:
        summary = build_clone_barrier_metrics(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clone_section_rows']} (clone, section) row(s), "
        f"{summary['n_distinct_clones']} distinct clone(s). Significant (p<0.05) barrier enrichment: "
        f"fibroblast={summary['n_significant_fibroblast_barrier']}, vascular={summary['n_significant_vascular_barrier']}, "
        f"suppressive_myeloid={summary['n_significant_suppressive_myeloid_barrier']}. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
