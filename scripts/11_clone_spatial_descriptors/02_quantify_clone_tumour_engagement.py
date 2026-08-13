#!/usr/bin/env python3
"""
`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`

Computes the central clone-tumour-engagement measures (malignant-cell
adjacency, opportunity-normalised engagement ratio, penetration,
interface localisation), cross-checked against Quality Control's
segmentation-artefact controls (spillover risk, resegmentation
concordance) -- see
src/xenium_tcr_ecology/clone_ecology/tumour_engagement.py's module
docstring.

Primary output: data/derived/clone_tumour_engagement.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.tumour_engagement import build_clone_tumour_engagement


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
        script_name="02_quantify_clone_tumour_engagement",
        project_root=project_root,
        phase="11_clone_spatial_descriptors",
    )

    try:
        summary = build_clone_tumour_engagement(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clone_section_rows']} (clone, section) row(s), "
        f"{summary['n_distinct_clones']} distinct clone(s), "
        f"{summary['n_rows_with_resegmentation_check']} with a resegmentation check. "
        f"Mean engagement ratio {summary['mean_engagement_ratio']:.3f}. "
        f"{summary['n_clone_sections_penetrating_tumour']} clone-section(s) penetrate the tumour mask. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
