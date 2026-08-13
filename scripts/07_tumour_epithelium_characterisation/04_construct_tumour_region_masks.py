#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`

Forms spatially coherent tumour regions from `07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`'s continuous
malignant-cell scores (thresholded at malignancy_score > 0, spatially
majority-vote smoothed over k=10 neighbours, then grouped into connected
components) and removes isolated false positives (regions smaller than
MIN_REGION_SIZE_CELLS). See
src/xenium_tcr_ecology/tumour/region_masks.py's module docstring for the
full method and threshold rationale.

Primary output: data/derived/tumour_masks/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.region_masks import build_tumour_region_masks


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "07_tumour_epithelium_characterisation",
        script_name="04_construct_tumour_region_masks",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_tumour_region_masks(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_sections']} section(s), {summary['n_cells_total']} cell(s). "
        f"{summary['n_in_tumour_region_total']} ({summary['fraction_in_tumour_region']*100:.2f}%) in a "
        f"tumour region across {summary['n_regions_total']} region(s); "
        f"{summary['n_isolated_removed_total']} isolated false positive(s) removed. "
        f"Wrote {summary['output_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
