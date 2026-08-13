#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`

Builds radius (15/30/50um), k-nearest-neighbour (k=6/10/15), Delaunay and
polygon-boundary-contact graphs for every section -- candidate
parameter values for `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` to later calibrate/select from, not yet a
single final choice. See
src/xenium_tcr_ecology/graphs/candidate_graphs.py's module docstring for
the full method and parameter grounding.

Primary output: data/graphs/candidate/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.candidate_graphs import build_all_candidate_graphs


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
        script_name="00_generate_candidate_spatial_graphs",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_all_candidate_graphs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_sections']} section(s), {summary['n_graph_types_per_section']} graph type(s) each "
        f"({summary['n_graphs_total']} graphs total). Median mean degree by type: "
        f"{summary['median_mean_degree_by_graph_type']}. Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
