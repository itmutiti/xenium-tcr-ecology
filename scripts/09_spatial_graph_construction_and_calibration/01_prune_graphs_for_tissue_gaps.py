#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`

Removes Delaunay/radius edges that bridge tissue folds, holes or gaps
between separated tissue fragments, using a per-section MAD-based
threshold derived from that section's Delaunay edge-length
distribution (this project's established outlier convention, Iglewicz &
Hoaglin 1993, already used identically in qc/spatial_artifacts.py). See
src/xenium_tcr_ecology/graphs/graph_pruning.py's module docstring.

Primary output: data/graphs/pruned/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.graph_pruning import prune_all_graphs


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
        script_name="01_prune_graphs_for_tissue_gaps",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = prune_all_graphs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_sections']} section(s). Median max edge length: "
        f"{summary['median_max_edge_length_um']}um. {summary['total_edges_pruned']} edge(s) pruned total. "
        f"Fraction pruned by type: {summary['fraction_edges_pruned_by_graph_type']}. "
        f"Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
