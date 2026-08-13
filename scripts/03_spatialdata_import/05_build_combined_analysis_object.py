#!/usr/bin/env python3
"""
`03_spatialdata_import/05_build_combined_analysis_object.py`

Concatenates all per-section AnnData objects (`03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` output, metadata
already attached) into one combined object with globally unique cell IDs
(section_id-prefixed) and the patient/run hierarchy already present in obs.
Gene panels are not identical across sections (patient-specific/batch-
specific CDR3 probes) -- uses an outer join and records per-section panel
membership separately so a 0-fill is never confused with "not measured".

Primary output: data/objects/hnscc_xenium_combined.h5ad;
                 results/tables/03_spatialdata_import/gene_panel_membership.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.combine_sections import build_combined_object


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    anndata_root = project_root / "data" / "objects" / "anndata"
    output_path = project_root / "data" / "objects" / "hnscc_xenium_combined.h5ad"
    panel_membership_path = (
        project_root
        / "results"
        / "tables"
        / "03_spatialdata_import"
        / "gene_panel_membership.parquet"
    )
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="05_build_combined_analysis_object",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    try:
        summary = build_combined_object(anndata_root, output_path, panel_membership_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   Combined {summary['n_sections']} section(s), {summary['n_patients']} patient(s): "
        f"{summary['n_cells_total']} cells x {summary['n_genes_union']} genes (union), "
        f"{summary['n_genes_core']} core genes shared by every section -> {output_path}"
    )
    print(f"[INFO] Per-section panel membership: {panel_membership_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
