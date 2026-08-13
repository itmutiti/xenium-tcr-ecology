#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`

Cross-checks the transcriptional malignancy call against morphology. No
pathologist tumour/normal region annotation exists in this dataset
(only raw DAPI/boundary-stain morphology images are
available), so a quantitative sensitivity/specificity concordance cannot
be computed here; see
src/xenium_tcr_ecology/tumour/morphology_concordance.py's module docstring
 for the two non-fabricated
partial checks this produces instead (spatial autocorrelation of the
malignancy score, and DAPI-image overlay panels for qualitative human
review), and why.

Primary output: reports/tumour/morphology_concordance.pdf
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.morphology_concordance import build_morphology_concordance_report


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
        script_name="02_cross_validate_against_morphology",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_morphology_concordance_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Median Moran's I of malignancy_probability across "
        f"{summary['n_sections_with_autocorrelation']} section(s): {summary['median_morans_i']:.4f} "
        f"(range {summary['min_morans_i']:.4f}-{summary['max_morans_i']:.4f}). "
        f"{summary['n_panels_rendered']} overlay panel(s) rendered ({summary['n_panels_failed']} failed). "
        f"Wrote {summary['autocorr_path']}, {summary['output_path']}. "
        f"STATUS: {summary['concordance_status']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
