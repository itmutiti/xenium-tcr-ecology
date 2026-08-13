#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`

Creates the prespecified primary graph (`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`'s calibrated
radius=30um, gap-pruned, upgraded to a distance-weighted graph) with
node metadata (`06_cell_type_annotation/06_integrate_annotation_evidence.py`'s final_lineage/final_substate/confidence)
and edge metadata (Euclidean distance in um), with patient-separated
connected components. See
src/xenium_tcr_ecology/graphs/primary_graph.py's module docstring for the radius-vs-kNN judgment call.

Primary output: data/graphs/primary_graphs/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.primary_graph import build_primary_graphs


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
        script_name="03_construct_primary_cell_graph",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_primary_graphs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Primary graph type: {summary['graph_type']}. {summary['n_sections']} section(s), "
        f"{summary['n_patients']} patient(s), {summary['n_nodes_total']} node(s), "
        f"{summary['n_edges_total']} edge(s). {summary['n_cross_patient_edges_found']} cross-patient "
        f"edge(s) found (must be 0). Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
