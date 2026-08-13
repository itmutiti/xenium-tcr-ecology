#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/04_construct_tumour_tcell_bipartite_graph.py`

Represents direct (touching cell-boundary polygons, `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`'s
boundary_contact graph) and near-direct (`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated 30um
radius graph) malignant-cell/T-cell relationships for contact-focused
analyses, extracted as cross-type subgraphs of the already-built,
already-validated `09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`, `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` graphs. See
src/xenium_tcr_ecology/graphs/tumour_tcell_bipartite.py's module
docstring for the malignant-cell population definition (`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`'s
spatially-validated in_tumour_region flag).

Primary output: data/graphs/tumour_tcell/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.tumour_tcell_bipartite import build_tumour_tcell_bipartite_graphs


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
        script_name="04_construct_tumour_tcell_bipartite_graph",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_tumour_tcell_bipartite_graphs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_sections']} section(s), {summary['n_malignant_cells_total']} malignant cell(s), "
        f"{summary['n_tcells_total']} T cell(s). Direct contact: {summary['n_direct_contact_edges_total']} edge(s), "
        f"{summary['n_tcells_with_any_direct_contact']} T cell(s) in contact. Near-direct: "
        f"{summary['n_near_direct_contact_edges_total']} edge(s), "
        f"{summary['n_tcells_with_any_near_direct_contact']} T cell(s) in contact. Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
