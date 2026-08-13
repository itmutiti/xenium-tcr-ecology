#!/usr/bin/env python3
"""
`05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`

Builds layers["counts"] (raw, preserved explicitly), layers["lognorm"]
(total-count normalisation + log1p), layers["pearson_residuals"] (analytic
Pearson residual variance stabilisation), and layers["detected"] (binary
detection) on `04_quality_control/07_apply_qc_filters_with_audit_trail.py`'s QC-filtered object. Both normalisation methods
compute per-cell exposure from biological_gene features only (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s
classification), not the full panel -- see
src/xenium_tcr_ecology/preprocess/count_layers.py's module docstring for
the composition-bias finding motivating this.

Primary output: data/objects/analysis_ready.h5ad
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.preprocess.count_layers import build_count_layers_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "05_preprocessing_and_normalisation",
        script_name="01_construct_analysis_count_layers",
        project_root=project_root,
        phase="05_preprocessing_and_normalisation",
    )

    try:
        summary = build_count_layers_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells x {summary['n_genes']} genes "
        f"({summary['n_exposure_genes']} exposure genes). "
        f"Layers: {summary['layers']}. Target sum: {summary['lognorm_target_sum']:.1f}, "
        f"median size factor: {summary['median_size_factor']:.3f}. "
        f"Wrote data/objects/analysis_ready.h5ad"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
