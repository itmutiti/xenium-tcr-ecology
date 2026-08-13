#!/usr/bin/env python3
"""
`14_spatial_interactions_and_barriers/04_analyse_barrier_pathways.py`

Evaluates checkpoint/chemokine/interferon/antigen-presentation programme
activity specifically at fibroblast/suppressive-myeloid barrier
interface cells, contrasting excluded (below-median engagement) vs
engaged (above-median engagement) T-cell clones. Exploratory -- follows
directly from `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s prespecified finding that
suppressive-myeloid barrier fraction predicts lower clone-tumour
engagement. See src/xenium_tcr_ecology/interactions/barrier_pathways.py's
module docstring.

Primary output: data/derived/barrier_pathways.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.interactions.barrier_pathways import build_barrier_pathways


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "14_spatial_interactions_and_barriers",
        script_name="04_analyse_barrier_pathways",
        project_root=project_root,
        phase="14_spatial_interactions_and_barriers",
    )

    try:
        summary = build_barrier_pathways(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_classified_clone_sections']} classified clone-section(s), "
        f"{summary['n_unique_interface_cells']} unique interface cell(s) "
        f"({summary['n_interface_cell_clone_rows']} cell-clone rows before dedup). "
        f"{summary['n_significant_bh']}/{summary['n_programs_tested']} programme(s) significant (BH q<0.05). "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
