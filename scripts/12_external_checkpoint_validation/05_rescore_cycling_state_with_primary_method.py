#!/usr/bin/env python3
"""
`12_external_checkpoint_validation/05_rescore_cycling_state_with_primary_method.py`

Re-scores GSE103322's T cells for the Cycling state using the pipeline's
primary scoring method (`scanpy.tl.score_genes`) instead of the
simplified z-score proxy used in `00_project_provisional_signatures_
to_bulk_reference.py`, to test whether the Cycling-state discrepancy
between this cohort (40.3%) and GSE103322 (14.4%)
reflects a cohort-composition difference or a scoring-method artefact.
Marker gene sets, seed, and the state-classification rule
(`xenium_tcr_ecology.external_checkpoint.bulk_reference.
assign_t_cell_substate`) are unchanged from the original comparison; only
the scoring method and gene pool (full transcriptome, 23,686 genes)
differ. Added after the pipeline was initially complete; see
`docs/analysis_amendments.md`.

Primary output: data/derived/cycling_rescore_comparison.parquet,
                 reports/external_checkpoint/cycling_rescore_comparison.pdf
"""

from __future__ import annotations

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from xenium_tcr_ecology.annotation.t_cell_substates import TREG_MARKERS
from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.external_checkpoint.bulk_reference import (
    assign_t_cell_substate,
    compare_state_proportions,
    parse_gse103322_full_matrix_t_cells,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.infra.seeding import get_annotation_seed
from xenium_tcr_ecology.preprocess.program_scores import PROGRAM_GENE_SETS

MODULE_SCORE_PROGRAMS = ["cytotoxicity", "exhaustion", "proliferation"]


def rescore_with_scanpy(expr_df: pd.DataFrame) -> pd.DataFrame:
    """Re-scores GSE103322 T cells using scanpy.tl.score_genes (the
    primary pipeline's scoring method) in place of the original
    simplified z-score proxy. Returns a DataFrame of per-cell scores,
    one column per program plus `treg_score`."""
    gene_cols = [c for c in expr_df.columns if c != "non-cancer cell type"]
    adata = ad.AnnData(
        X=expr_df[gene_cols].to_numpy(dtype=np.float32),
        obs=expr_df[["non-cancer cell type"]].copy(),
    )
    adata.var_names = gene_cols
    adata.layers["log_norm"] = adata.X.copy()

    gene_pool = gene_cols  # full-transcriptome pool: every gene GSE103322 measures

    seed = get_annotation_seed()
    scores = pd.DataFrame(index=adata.obs_names)
    for program in MODULE_SCORE_PROGRAMS:
        genes = [g for g in PROGRAM_GENE_SETS[program] if g in adata.var_names]
        sc.tl.score_genes(
            adata,
            gene_list=genes,
            gene_pool=gene_pool,
            layer="log_norm",
            score_name=f"{program}_score",
            random_state=seed,
        )
        scores[f"{program}_score"] = adata.obs[f"{program}_score"]

    treg_genes = [g for g in TREG_MARKERS if g in adata.var_names]
    sc.tl.score_genes(
        adata,
        gene_list=treg_genes,
        gene_pool=gene_pool,
        layer="log_norm",
        score_name="treg_score",
        random_state=seed,
    )
    scores["treg_score"] = adata.obs["treg_score"]

    scores["cd4_expr"] = expr_df["CD4"].values if "CD4" in expr_df.columns else np.nan
    scores["cd8a_expr"] = expr_df["CD8A"].values if "CD8A" in expr_df.columns else np.nan
    return scores


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "12_external_checkpoint_validation",
        script_name="05_rescore_cycling_state_with_primary_method",
        project_root=project_root,
        phase="12_external_checkpoint_validation",
    )

    gse103322_path = (
        project_root / "data" / "external" / "GSE103322" / "GSE103322_HNSCC_all_data.txt.gz"
    )
    t_cell_states_path = project_root / "data" / "derived" / "t_cell_states.parquet"
    original_comparison_path = (
        project_root / "data" / "derived" / "bulk_reference_state_comparison.parquet"
    )

    for p in (gse103322_path, t_cell_states_path, original_comparison_path):
        if not p.is_file():
            print(
                f"[ERROR] '{p}' not found. Run `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py` first.",
                file=sys.stderr,
            )
            return 1

    try:
        expr_df = parse_gse103322_full_matrix_t_cells(gse103322_path)
        scanpy_scores = rescore_with_scanpy(expr_df)
        scanpy_states = pd.Series(
            assign_t_cell_substate(
                scanpy_scores["treg_score"].to_numpy(),
                scanpy_scores["proliferation_score"].to_numpy(),
                scanpy_scores["exhaustion_score"].to_numpy(),
                scanpy_scores["cytotoxicity_score"].to_numpy(),
                scanpy_scores["cd4_expr"].to_numpy(),
                scanpy_scores["cd8a_expr"].to_numpy(),
            ),
            index=expr_df.index,
        )

        project_states = pd.read_parquet(t_cell_states_path)["t_cell_state"]
        scanpy_comparison = compare_state_proportions(scanpy_states, project_states)
        scanpy_comparison = scanpy_comparison.rename(
            columns={"reference_fraction": "reference_fraction_scanpy_rescore"}
        )

        original_comparison = pd.read_parquet(original_comparison_path).rename(
            columns={"reference_fraction": "reference_fraction_original_zscore_proxy"}
        )

        comparison = original_comparison.merge(
            scanpy_comparison[["state", "reference_fraction_scanpy_rescore"]],
            on="state",
            how="outer",
        )
        comparison = comparison[
            [
                "state",
                "reference_fraction_original_zscore_proxy",
                "reference_fraction_scanpy_rescore",
                "project_fraction",
            ]
        ]
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    output_path = project_root / "data" / "derived" / "cycling_rescore_comparison.parquet"
    comparison.to_parquet(output_path)

    report_dir = project_root / "reports" / "external_checkpoint"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / "cycling_rescore_comparison.pdf"

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(comparison))
    width = 0.25
    ax.bar(
        x - width,
        comparison["reference_fraction_original_zscore_proxy"],
        width,
        label="reference, original z-score proxy",
    )
    ax.bar(
        x,
        comparison["reference_fraction_scanpy_rescore"],
        width,
        label="reference, scanpy.tl.score_genes re-score",
    )
    ax.bar(x + width, comparison["project_fraction"], width, label="this project (primary cohort)")
    ax.set_xticks(x)
    ax.set_xticklabels(comparison["state"], rotation=30, ha="right")
    ax.set_ylabel("Fraction of T cells")
    ax.set_title(
        "T-cell state proportions: does the scoring method explain the Cycling discrepancy?"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(report_path)
    plt.close(fig)

    cycling_row = comparison[comparison["state"] == "Cycling"].iloc[0]
    logger.log_event(
        n_t_cells=len(expr_df),
        cycling_fraction_original_proxy=float(
            cycling_row["reference_fraction_original_zscore_proxy"]
        ),
        cycling_fraction_scanpy_rescore=float(cycling_row["reference_fraction_scanpy_rescore"]),
        cycling_fraction_project=float(cycling_row["project_fraction"]),
        output_path=str(output_path),
    )
    logger.write(status="ok")

    print(f"[OK]   Re-scored {len(expr_df)} GSE103322 T cells with scanpy.tl.score_genes.")
    print(
        f"[INFO] Cycling fraction: original z-score proxy={cycling_row['reference_fraction_original_zscore_proxy']:.3f}, scanpy re-score={cycling_row['reference_fraction_scanpy_rescore']:.3f}, this project's primary cohort={cycling_row['project_fraction']:.3f}"
    )
    print(f"[OK]   Wrote {output_path}, {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
