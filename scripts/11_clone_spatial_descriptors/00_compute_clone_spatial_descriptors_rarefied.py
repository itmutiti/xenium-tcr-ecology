#!/usr/bin/env python3
"""
`11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`

Computes five rarefaction-normalised clone spatial descriptors
(dispersion, clustering, tumour distance, border enrichment, domain
occupancy) per (clone, section), restricted to the primary HNSCC cohort
-- see src/xenium_tcr_ecology/clone_ecology/spatial_descriptors.py's
module docstring.

Primary output: data/derived/clone_spatial_descriptors.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.spatial_descriptors import build_clone_spatial_descriptors


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "11_clone_spatial_descriptors",
        script_name="00_compute_clone_spatial_descriptors_rarefied",
        project_root=project_root,
        phase="11_clone_spatial_descriptors",
    )

    try:
        summary = build_clone_spatial_descriptors(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clone_section_rows']} (clone, section) row(s), "
        f"{summary['n_distinct_clones']} distinct clone(s), {summary['n_sections']} section(s). "
        f"{summary['n_rarefied_dispersion']} row(s) with rarefied descriptors (>=5 cells). "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
