#!/usr/bin/env python3
"""
`03_spatialdata_import/03_create_anndata_expression_objects.py`

Extracts each section's expression table from its SpatialData Zarr store
into a standalone AnnData .h5ad, preserving the spatial link (obsm['spatial'],
cell_id) so downstream scripts that only need expression data don't have to
load the full spatial store.

Primary output: data/objects/anndata/<section>.h5ad
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.anndata_export import export_section_anndata


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_root = project_root / "data" / "objects" / "anndata"
    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        print(
            f"[ERROR] No .zarr stores found under '{spatialdata_root}'. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first.",
            file=sys.stderr,
        )
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="03_create_anndata_expression_objects",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    exported, failed = 0, []
    for zarr_path in zarr_paths:
        output_path = output_root / f"{zarr_path.stem}.h5ad"
        try:
            summary = export_section_anndata(zarr_path, output_path)
        except PipelineError as exc:
            print(f"[ERROR] {zarr_path.stem}: {exc}", file=sys.stderr)
            logger.log_error(f"{zarr_path.stem}: {exc}")
            failed.append(zarr_path.stem)
            continue
        print(
            f"[OK]   {summary['section_id']}: {summary['n_cells']} cells x {summary['n_genes']} genes -> {output_path}"
        )
        exported += 1

    status = "ok" if not failed else "failed"
    logger.log_event(exported=exported, failed=failed)
    logger.write(status=status)
    print(f"\n[SUMMARY] {exported} exported, {len(failed)} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
