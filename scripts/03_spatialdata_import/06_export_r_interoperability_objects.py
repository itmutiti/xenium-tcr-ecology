#!/usr/bin/env python3
"""
`03_spatialdata_import/06_export_r_interoperability_objects.py`

Exports each section's AnnData as a standard 10x Matrix Market triplet
(matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz, readable natively by
Seurat's Read10X()) plus obs/var parquet carrying spatial coordinates and
all clinical/technical metadata -- deliberately not using a direct h5ad ->
Seurat conversion tool (see module docstring in
xenium_tcr_ecology.io.r_export for why).

Primary output: data/objects/r_exports/<section>/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.r_export import export_section_r_interop


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    anndata_root = project_root / "data" / "objects" / "anndata"
    output_root = project_root / "data" / "objects" / "r_exports"
    h5ad_paths = sorted(anndata_root.glob("*.h5ad"))
    if not h5ad_paths:
        print(
            f"[ERROR] No .h5ad files found under '{anndata_root}'. Run `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` first.",
            file=sys.stderr,
        )
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="06_export_r_interoperability_objects",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    exported, failed = 0, []
    for h5ad_path in h5ad_paths:
        section_id = h5ad_path.stem
        try:
            summary = export_section_r_interop(h5ad_path, output_root / section_id)
        except PipelineError as exc:
            print(f"[ERROR] {section_id}: {exc}", file=sys.stderr)
            logger.log_error(f"{section_id}: {exc}")
            failed.append(section_id)
            continue
        print(
            f"[OK]   {section_id}: {summary['n_cells']} cells x {summary['n_genes']} genes -> {summary['output_dir']}"
        )
        exported += 1

    status = "ok" if not failed else "failed"
    logger.log_event(exported=exported, failed=failed)
    logger.write(status=status)
    print(f"\n[SUMMARY] {exported} exported, {len(failed)} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
