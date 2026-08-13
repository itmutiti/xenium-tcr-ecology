#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`

Combines tumour keratins, HPV genes, stress/EMT programs and a
patient-clonality proxy (substituting for external reference mapping --
see src/xenium_tcr_ecology/tumour/malignancy_scoring.py's module docstring
 for why: this project's only
whole-transcriptome external reference, GSE287301, is T-cell-only) to
estimate a continuous per-cell malignancy probability within `07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py`'s
epithelial subset.

Primary output: data/derived/malignancy_scores.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.malignancy_scoring import build_malignancy_score_report


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
        script_name="01_score_malignancy_and_normal_epithelium",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_malignancy_score_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']} epithelial cell(s) scored. "
        f"{summary['n_cells_hpv_scored']} ({summary['fraction_cells_hpv_scored']*100:.1f}%) "
        f"had an HPV-probe panel available. Mean malignancy score: {summary['mean_malignancy_score']:.4f}. "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
