#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`

Creates tumour-border geometries (a buffered-union polygon per section
from `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`'s in-region cells) and signed distance-to-boundary values
for ALL cells in the primary analysis matrix, not only epithelial cells --
see src/xenium_tcr_ecology/tumour/tumour_boundaries.py's module docstring
for the full method and the signed-distance sign convention.

Primary output: data/derived/tumour_boundaries.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.tumour_boundaries import build_tumour_boundaries


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
        script_name="05_extract_tumour_boundaries",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_tumour_boundaries(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} cell(s) across {summary['n_sections']} section(s) "
        f"({summary['n_sections_with_tumour_mask']} with a tumour mask). "
        f"{summary['n_cells_inside_tumour_region_total']} "
        f"({summary['fraction_cells_inside_tumour_region']*100:.2f}%) inside a tumour region. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
