#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`

Selects one calibrated radius and one calibrated k from `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`'s
candidates, using connected-component analysis on the `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`
pruned graphs across all 18 sections -- see
src/xenium_tcr_ecology/graphs/graph_calibration.py's module docstring
rule.

Primary output: config/graph_parameters.yaml
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.graph_calibration import build_graph_parameter_calibration


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
        script_name="02_calibrate_graph_parameters",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_graph_parameter_calibration(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Calibrated radius: {summary['calibrated_radius_um']}um, calibrated k: {summary['calibrated_knn']}. "
        f"Median largest-component fraction by candidate: {summary['median_largest_component_fraction_by_candidate']}. "
        f"{summary['n_sections_with_genuine_multi_fragment_structure']} section(s) show genuine multi-fragment "
        f"tissue structure. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
