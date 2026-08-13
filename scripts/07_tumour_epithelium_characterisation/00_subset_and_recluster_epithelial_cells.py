#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py`

Characterises epithelial heterogeneity separately from immune/stromal
structure: subsets `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s `Epithelial_Tumour`-lineage cells and
reclusters them at both joint and within-patient scope, reusing Phase
6.01's clustering functions. See
src/xenium_tcr_ecology/tumour/epithelial_subset.py's module docstring for
the full method and the expected (literature-grounded) patient-dominance
finding this milestone measures and records.

Primary output: data/objects/epithelial_subset.h5ad
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.epithelial_subset import build_epithelial_subset_report


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
        script_name="00_subset_and_recluster_epithelial_cells",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = build_epithelial_subset_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_epithelial_cells']} epithelial cell(s) across "
        f"{summary['n_patients_represented']} patient(s). "
        f"Joint-cluster patient dominance (res={summary['diagnostic_resolution']}): "
        f"{summary['joint_cluster_patient_dominance']:.4f} (unweighted), "
        f"{summary['joint_cluster_patient_dominance_cell_weighted']:.4f} (cell-weighted). "
        f"Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
