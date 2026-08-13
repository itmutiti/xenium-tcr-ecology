"""Export each section's AnnData as R-interoperable files (`03_spatialdata_import/06_export_r_interoperability_objects.py`).

Deliberately NOT using a direct h5ad->Seurat conversion tool (e.g.
SeuratDisk): those are version-fragile and less actively maintained than
either side of the conversion they bridge. Instead exports the standard 10x
Matrix Market triplet (matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz),
readable natively by Seurat's own Read10X() with zero extra R packages,
plus separate obs/var parquet files carrying everything Read10X does not
(spatial coordinates, patient/clinical metadata, gene metadata) so no
information is lost, only reshaped into a more universally-portable form.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from xenium_tcr_ecology.infra.exceptions import PipelineError


def export_section_r_interop(anndata_path: Path, output_dir: Path) -> dict:
    if not anndata_path.is_file():
        raise PipelineError(
            f"'{anndata_path}' not found. Run `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` first."
        )

    adata = ad.read_h5ad(anndata_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Matrix Market expects features (genes) as rows, cells as columns --
    # the transpose of AnnData's cells x genes convention.
    matrix = csr_matrix(adata.X).T
    with gzip.open(output_dir / "matrix.mtx.gz", "wb") as f:
        mmwrite(f, matrix)

    with gzip.open(output_dir / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(adata.obs_names) + "\n")

    with gzip.open(output_dir / "features.tsv.gz", "wt") as f:
        for gene in adata.var_names:
            f.write(f"{gene}\t{gene}\tGene Expression\n")

    obs = adata.obs.copy()
    if "spatial" in adata.obsm:
        obs["x_centroid"] = adata.obsm["spatial"][:, 0]
        obs["y_centroid"] = adata.obsm["spatial"][:, 1]
    obs.to_parquet(output_dir / "obs_metadata.parquet")
    adata.var.to_parquet(output_dir / "var_metadata.parquet")

    return {"n_cells": adata.n_obs, "n_genes": adata.n_vars, "output_dir": str(output_dir)}
