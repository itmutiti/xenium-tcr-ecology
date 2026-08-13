#!/usr/bin/env python3
"""
`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`

Computes clone-level proximity and opportunity-normalised enrichment
with dendritic cells, macrophages, and antigen-presentation
programme activity in a clone's spatial neighbourhood -- see
src/xenium_tcr_ecology/clone_ecology/apc_support.py's module docstring
.

Primary output: data/derived/clone_apc_support.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.apc_support import build_clone_apc_support


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
        script_name="03_quantify_clone_apc_support",
        project_root=project_root,
        phase="11_clone_spatial_descriptors",
    )

    try:
        summary = build_clone_apc_support(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clone_section_rows']} (clone, section) row(s), "
        f"{summary['n_distinct_clones']} distinct clone(s). Mean DC/Macrophage engagement ratios: "
        f"DC={summary['mean_dc_engagement_ratio']:.3f}, Macrophage={summary['mean_macrophage_engagement_ratio']:.3f}. "
        f"Mean antigen-presentation score excess={summary['mean_antigen_presentation_score_excess']:.4f}. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
