#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/05_construct_clone_induced_subgraphs.py`

Extracts all cells belonging to each `08_tcr_clonal_analysis/08_generate_tcr_release_report.py` high-confidence clone
plus 3 concentric microenvironment shells (BFS graph-distance rings over
`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated primary graph). See
src/xenium_tcr_ecology/graphs/clone_subgraphs.py's module docstring.

Primary output: data/graphs/clones/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.clone_subgraphs import build_clone_subgraphs


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
        script_name="05_construct_clone_induced_subgraphs",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_clone_subgraphs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clones']} clone(s), {summary['n_clone_section_instances']} clone-section instance(s). "
        f"Median clone size: {summary['median_n_clone_cells']}, median shell-1 size: "
        f"{summary['median_n_shell1_cells']}. Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
