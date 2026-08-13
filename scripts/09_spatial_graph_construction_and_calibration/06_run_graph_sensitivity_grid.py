#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/06_run_graph_sensitivity_grid.py`

Repeats two core metrics (graph connectivity, tumour-T-cell contact rate)
over all six `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py` candidate radii/k values, and checks whether the
tumour-T-cell contact rate's per-section conclusions are robust to the
exact radius chosen (Spearman rank correlation against `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`, `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s
calibrated 30um choice). See
src/xenium_tcr_ecology/graphs/sensitivity_grid.py's module docstring.

Primary output: reports/graphs/sensitivity_grid.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.sensitivity_grid import build_sensitivity_grid


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root
        / "results"
        / "logs"
        / "09_spatial_graph_construction_and_calibration",
        script_name="06_run_graph_sensitivity_grid",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_sensitivity_grid(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_sections']} section(s), {summary['n_candidate_graph_types']} graph type(s), "
        f"{summary['n_rows']} row(s). Median fraction T cells in contact by graph: "
        f"{summary['median_fraction_tcells_in_contact_by_graph']}. Spearman rho vs. calibrated radius: "
        f"{summary['spearman_rho_vs_calibrated_radius']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
