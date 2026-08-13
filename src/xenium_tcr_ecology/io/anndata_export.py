"""Extract each section's expression table from its SpatialData store into a
standalone AnnData .h5ad (`03_spatialdata_import/03_create_anndata_expression_objects.py`), preserving the link back to spatial
elements via obsm['spatial'] and a retained cell_id / section_id column
(rather than exporting counts alone and losing spatial context)."""

from __future__ import annotations

from pathlib import Path

import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError


def export_section_anndata(zarr_path: Path, output_path: Path) -> dict:
    if not zarr_path.exists():
        raise PipelineError(
            f"'{zarr_path}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    sdata = sd.read_zarr(zarr_path)
    adata = sdata["table"].copy()

    if "spatial" not in adata.obsm:
        raise PipelineError(
            f"'{zarr_path}': table has no obsm['spatial'] -- spatial link would be lost."
        )
    if "cell_id" not in adata.obs.columns:
        raise PipelineError(
            f"'{zarr_path}': table has no obs['cell_id'] -- spatial element link would be lost."
        )

    section_id = zarr_path.stem
    adata.obs["section_id"] = section_id

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return {"section_id": section_id, "n_cells": adata.n_obs, "n_genes": adata.n_vars}
