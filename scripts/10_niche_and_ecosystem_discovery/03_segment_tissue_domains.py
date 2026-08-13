#!/usr/bin/env python3
"""
`10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`

Forms spatially contiguous tissue domains from `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s per-cell
archetype labels via majority-vote smoothing over `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated
30um primary graph, followed by same-label connected components -- see
src/xenium_tcr_ecology/niches/tissue_domains.py's module docstring.

Primary output: data/derived/tissue_domains.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.niches.tissue_domains import build_tissue_domains


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "10_niche_and_ecosystem_discovery",
        script_name="03_segment_tissue_domains",
        project_root=project_root,
        phase="10_niche_and_ecosystem_discovery",
    )

    try:
        summary = build_tissue_domains(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} cell(s) ({summary['n_cells_excluded_zero_degree']} zero-degree excluded), "
        f"{summary['n_cells_relabelled_by_smoothing']} relabelled by majority-vote smoothing, "
        f"{summary['n_domains']} domain(s) ({summary['n_single_cell_domains']} single-cell). "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
