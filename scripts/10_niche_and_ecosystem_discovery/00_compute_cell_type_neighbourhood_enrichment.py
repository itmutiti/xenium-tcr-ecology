#!/usr/bin/env python3
"""
`10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py`

Tests pairwise cell-type adjacency using within-section constrained
permutations (calibrated in `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`), via squidpy's nhood_enrichment
on `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated primary graph. See
src/xenium_tcr_ecology/niches/neighbourhood_enrichment.py's module
docstring for the full method.

Primary output: data/derived/neighbourhood_enrichment.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.niches.neighbourhood_enrichment import build_neighbourhood_enrichment


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
        script_name="00_compute_cell_type_neighbourhood_enrichment",
        project_root=project_root,
        phase="10_niche_and_ecosystem_discovery",
    )

    try:
        summary = build_neighbourhood_enrichment(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_sections']} section(s), {summary['n_lineage_pairs']} lineage pair(s). "
        f"Top enriched: {summary['top_enriched_pairs']}. Top depleted: {summary['top_depleted_pairs']}. "
        f"Wrote {summary['output_path']}, {summary['summary_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
