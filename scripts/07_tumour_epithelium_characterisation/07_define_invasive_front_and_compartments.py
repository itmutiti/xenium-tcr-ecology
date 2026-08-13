#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/07_define_invasive_front_and_compartments.py`

Labels tumour core, inner margin, outer margin and distal stroma from
`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s signed distance-to-boundary, using a predeclared band width
grounded directly in this dataset's achievable interior depth (not
literature values that do not fit this dataset's boundary geometry -- see
src/xenium_tcr_ecology/tumour/spatial_compartments.py's module docstring
), plus two sensitivity band
widths per the blueprint's explicit requirement.

Primary output: data/derived/spatial_compartments.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.spatial_compartments import build_spatial_compartments


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
        script_name="07_define_invasive_front_and_compartments",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_spatial_compartments(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} cell(s), {summary['n_cells_with_compartment']} with a compartment "
        f"assignment ({summary['n_cells_no_tumour_mask']} in sections with no tumour mask). "
        f"Primary band width {summary['primary_band_width_um']}um: {summary['primary_compartment_counts']}. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
