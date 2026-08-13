"""T-cell substate input preparation (`06_cell_type_annotation/04_resolve_t_cell_substates.R` helper).

Computes the marker evidence `04_resolve_t_cell_substates.R` needs to
assign discrete T-cell substates (CD4, CD8, Treg, Cycling, Cytotoxic,
Exhausted, Ambiguous), and writes it to a plain parquet the R script can
read directly (unlike analysis_ready.h5ad, which R cannot read -- the same
"Python computes, R decides/reports" split already established in Phase
5.02).

Reuses `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s `cytotoxicity_score`, `exhaustion_score` and
`proliferation_score` (program_scores.parquet) directly rather than
recomputing them -- the same marker sets would otherwise be duplicated.
Only Treg scoring and direct CD4/CD8A expression are new here: no existing
`05_preprocessing_and_normalisation/03_calculate_program_scores.py` program covers Treg identity, and CD4-vs-CD8 lineage is a direct
marker comparison, not a multi-gene program score.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

# Confirmed present in the panel; IKZF2 (also a canonical Treg marker)
# is not in this panel and is therefore not included.
TREG_MARKERS = ["FOXP3", "IL2RA", "CTLA4"]


def compute_treg_score(
    adata: ad.AnnData, layer: str, gene_pool: list[str], rng_seed: int = get_annotation_seed()
) -> pd.Series:
    present = [g for g in TREG_MARKERS if g in adata.var_names]
    if len(present) < 2:
        raise PipelineError(f"Only {len(present)} Treg marker(s) present in adata.var_names.")
    sc.tl.score_genes(
        adata,
        gene_list=present,
        gene_pool=gene_pool,
        layer=layer,
        score_name="treg_score",
        random_state=rng_seed,
    )
    return adata.obs["treg_score"]


def build_t_cell_substate_inputs(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    program_scores_path = project_root / "data" / "derived" / "program_scores.parquet"
    lineage_scores_path = project_root / "data" / "derived" / "lineage_scores.parquet"
    reference_labels_path = project_root / "data" / "derived" / "reference_labels.parquet"
    output_path = project_root / "data" / "derived" / "t_cell_substate_inputs.parquet"

    for p in (matrix_path, program_scores_path, lineage_scores_path, reference_labels_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    gene_pool = adata.var_names[adata.var["is_exposure_gene"]].tolist()

    treg_score = compute_treg_score(adata, layer, gene_pool)

    cd4_idx = adata.var_names.get_loc("CD4")
    cd8a_idx = adata.var_names.get_loc("CD8A")
    expr = adata.layers[layer]
    cd4_expr = expr[:, cd4_idx]
    cd8a_expr = expr[:, cd8a_idx]
    cd4_expr = (
        cd4_expr.toarray().ravel()
        if sparse.issparse(cd4_expr)
        else pd.Series(cd4_expr).to_numpy().ravel()
    )
    cd8a_expr = (
        cd8a_expr.toarray().ravel()
        if sparse.issparse(cd8a_expr)
        else pd.Series(cd8a_expr).to_numpy().ravel()
    )

    result = pd.DataFrame(
        {
            "treg_score": treg_score.to_numpy(),
            "cd4_expr": cd4_expr,
            "cd8a_expr": cd8a_expr,
        },
        index=adata.obs_names,
    )

    program_scores = pd.read_parquet(program_scores_path)
    result = result.join(
        program_scores[["cytotoxicity_score", "exhaustion_score", "proliferation_score"]]
    )

    lineage_scores = pd.read_parquet(lineage_scores_path)
    lineage_cols = [c for c in lineage_scores.columns if c.endswith("_lineage_score")]
    result["argmax_lineage"] = (
        lineage_scores[lineage_cols].idxmax(axis=1).str.replace("_lineage_score", "", regex=False)
    )

    reference_labels = pd.read_parquet(reference_labels_path)
    result = result.join(
        reference_labels[["predicted_state", "confidence"]].rename(
            columns={
                "predicted_state": "reference_predicted_state",
                "confidence": "reference_confidence",
            }
        )
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_cells": len(result),
        "n_t_or_nk_lineage_cells": int(result["argmax_lineage"].isin(["T_cell", "NK_cell"]).sum()),
    }
