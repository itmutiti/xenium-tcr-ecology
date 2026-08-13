#!/usr/bin/env python3
"""
`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`

Creates cell-centred neighbourhood composition vectors (fraction of each
of the 12 major lineages among spatial neighbours, own lineage
excluded) at Spatial Graph Construction and Calibration's three established candidate scales (radius 15/30/
50um, gap-pruned). See
src/xenium_tcr_ecology/niches/local_composition.py's module docstring.

Primary output: data/derived/local_compositions.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.niches.local_composition import build_local_compositions


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
        script_name="01_compute_local_neighbourhood_compositions",
        project_root=project_root,
        phase="10_niche_and_ecosystem_discovery",
    )

    try:
        summary = build_local_compositions(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} cell(s), {summary['n_lineages']} lineage(s), {summary['n_scales']} scale(s), "
        f"{summary['n_composition_columns']} composition column(s). Zero-degree cells by scale: "
        f"{summary['n_cells_zero_degree_by_scale']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
