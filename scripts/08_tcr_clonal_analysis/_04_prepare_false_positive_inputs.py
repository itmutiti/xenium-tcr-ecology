#!/usr/bin/env python3
"""
`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R` helper -- _04_prepare_false_positive_inputs.py

NOT its own numbered blueprint phase step: computes and exports the three
empirical negative controls (off-patient, non-T-cell, spatial
autocorrelation) that 04_estimate_false_positive_tcr_calls.R needs, since
that R script cannot read .h5ad directly (no R HDF5/AnnData reader
available). Invoked
by the R script via a subprocess call, not run standalone as a pipeline
phase.

Primary output: data/derived/tcr_false_positive_controls.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tcr.false_positive_estimation import prepare_false_positive_inputs


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
        script_name="_04_prepare_false_positive_inputs",
        project_root=project_root,
        phase="08_tcr_clonal_analysis",
    )

    try:
        summary = prepare_false_positive_inputs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_probes']} probe(s), {summary['n_probes_with_spatial_test']} with a spatial "
        f"autocorrelation test. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
