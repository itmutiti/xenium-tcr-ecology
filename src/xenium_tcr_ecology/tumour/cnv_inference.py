"""CNV-appendix input preparation (`07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R` helper).

`03_infer_cnv_appendix_only.R` cannot read `.h5ad` directly (no R HDF5/
AnnData reader available),
so this module -- invoked from R via `system2()`, matching the
`_02_compute_normalization_benchmark_metrics.py` precedent (`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`) --
exports what R needs as parquet:

1. **Epithelial cell-level expression** for the panel genes with resolved
   genomic coordinates (`07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R`'s own `gene_coordinates.py`; 395/399
   `biological_gene` features resolved, 4 dropped as unresolved
   rather than guessed).
2. **A per-patient reference baseline** -- mean lognorm expression per gene
   among that patient's own non-epithelial (immune/stromal) cells, not a
   cell-level export of the full reference population (723,533 cells
   across all patients; only the per-patient MEAN is actually needed
   downstream, so exporting it directly avoids moving a needless amount of
   data across the Python/R boundary). Patient-MATCHED reference cells
   are used, not a pooled cross-patient reference: this is standard CNV-
   inference practice (inferCNV/CopyKAT) and avoids introducing a
   cross-patient technical-batch confound into the very comparison this
   analysis exists to make.

Reference cells are non-epithelial (`06_cell_type_annotation/06_integrate_annotation_evidence.py`'s `final_lineage
!= "Epithelial_Tumour"`), not a "low malignancy_score" epithelial subset:
using the latter would be circular (the reference for validating the
malignancy score would itself be defined by that score).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.gene_coordinates import build_gene_coordinate_reference


def _to_dense_frame(adata: ad.AnnData, genes: list[str], layer: str) -> pd.DataFrame:
    X = adata[:, genes].layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    return pd.DataFrame(np.asarray(X), columns=genes, index=adata.obs_names)


def prepare_cnv_inputs(project_root: Path) -> dict:
    epithelial_subset_path = project_root / "data" / "objects" / "epithelial_subset.h5ad"
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    malignancy_scores_path = project_root / "data" / "derived" / "malignancy_scores.parquet"

    epithelial_expr_path = project_root / "data" / "derived" / "cnv_epithelial_expression.parquet"
    reference_baseline_path = project_root / "data" / "derived" / "cnv_reference_baseline.parquet"

    for p in (epithelial_subset_path, matrix_path, final_annotations_path, malignancy_scores_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    coord_summary = build_gene_coordinate_reference(project_root)
    coords = pd.read_csv(project_root / "references" / "gene_genomic_coordinates.tsv", sep="\t")

    sub = ad.read_h5ad(epithelial_subset_path)
    layer = sub.uns["primary_normalization_layer"]
    genes = [g for g in coords["gene"] if g in sub.var_names]
    if len(genes) == 0:
        raise PipelineError(
            "No coordinate-resolved genes found in epithelial_subset.h5ad var_names."
        )

    malignancy = pd.read_parquet(malignancy_scores_path)

    epi_expr = _to_dense_frame(sub, genes, layer)
    epi_expr["patient_id"] = sub.obs["patient_id"].to_numpy()
    epi_expr["malignancy_probability"] = malignancy.reindex(epi_expr.index)[
        "malignancy_probability"
    ].to_numpy()
    epi_expr.to_parquet(epithelial_expr_path)

    final_annotations = pd.read_parquet(final_annotations_path)
    reference_ids = final_annotations.index[
        final_annotations["final_lineage"] != "Epithelial_Tumour"
    ]

    primary = ad.read_h5ad(matrix_path)
    reference_ids = primary.obs_names.intersection(reference_ids)
    ref_adata = primary[reference_ids]
    ref_layer = primary.uns["primary_normalization_layer"]
    ref_genes = [g for g in genes if g in ref_adata.var_names]
    ref_expr = _to_dense_frame(ref_adata, ref_genes, ref_layer)
    ref_expr["patient_id"] = ref_adata.obs["patient_id"].to_numpy()
    reference_baseline = ref_expr.groupby("patient_id", observed=True)[ref_genes].mean()
    # Explicit reset_index(), not relying on the parquet writer to preserve
    # a named pandas index as a readable column across the R boundary --
    # confirmed this does not round-trip as a "patient_id" column
    # via arrow::read_parquet() in R (silently becomes an anonymous row
    # index instead), the same general class of Python/R interop hazard
    # already documented for booleans.
    reference_baseline = reference_baseline.reset_index()
    reference_baseline.to_parquet(reference_baseline_path)

    return {
        "n_genes_used": len(genes),
        "n_epithelial_cells": len(epi_expr),
        "n_reference_cells": len(ref_expr),
        "n_patients_with_reference_baseline": len(reference_baseline),
        "gene_coordinate_summary": coord_summary,
        "epithelial_expr_path": str(epithelial_expr_path),
        "reference_baseline_path": str(reference_baseline_path),
    }
