"""Transcriptional program transfer test (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`).

Directly investigates `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`'s finding: `cycling_fraction` -- the
dominant driver of `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s continuous ecological-structure
axis -- shows the largest state-abundance discrepancy between this
project's data and the independent GSE103322 reference (40.3% vs.
14.4%). This milestone asks a different, more fundamental question
than `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`'s abundance comparison: is each programme's
underlying gene-expression module itself an internally coherent,
co-expressed unit in both datasets, independent of any classification
threshold or platform-specific abundance difference?

**Coherence statistic:** mean pairwise Spearman correlation among a
programme's marker genes, within T cells, in each dataset separately.
Compared against an empirical null (mean pairwise correlation among
`N_PERMUTATIONS=200` random equal-sized gene sets drawn from the same
expressed-gene pool) -- a standard "is this gene set more coherently
co-expressed than chance" module-coherence test, not a comparison of
raw correlation magnitudes across platforms (which would not be a fair
comparison given Smart-seq2's far higher per-cell gene detection than a
623-gene targeted panel).

If a programme's genes are coherently co-expressed modules in both
datasets (both coherence values clear their own dataset's null), the
programme transfers as a biological unit -- any abundance discrepancy
found in `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py` then more likely reflects cohort
biology or a scoring-threshold difference, not a broken/non-transferable
gene signature. If a programme fails to clear its null in the reference
specifically, that is direct evidence the module itself does not
reproduce externally.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from xenium_tcr_ecology.annotation.t_cell_substates import TREG_MARKERS
from xenium_tcr_ecology.external_checkpoint.bulk_reference import (
    list_gse103322_genes,
    parse_gse103322_t_cells,
    required_gse103322_genes,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed
from xenium_tcr_ecology.preprocess.program_scores import PROGRAM_GENE_SETS

PROGRAMS_TESTED = ["cytotoxicity", "exhaustion", "proliferation"]
N_PERMUTATIONS = 200
RNG_SEED = get_default_seed()
N_BACKGROUND_GENES = 400


def filter_expressed_genes(
    expr_df: pd.DataFrame, genes: list[str], min_variance: float = 1e-6
) -> list[str]:
    """Restricts `genes` to those with non-degenerate variance in
    `expr_df`. A randomly sampled background gene that is never
    expressed (exactly constant, e.g. zero, across every cell) has an
    undefined correlation with anything -- including it in a null-model
    draw wastes a permutation and, in the unlucky case where an entire
    draw is degenerate, produces an all-NaN correlation matrix (on the
    first run: 15.5% of a 400-gene random background sample from
    GSE103322 had exactly zero variance among the 1,237 reference T
    cells)."""
    variances = expr_df[genes].var()
    return variances[variances >= min_variance].index.tolist()


def compute_mean_pairwise_correlation(expr_df: pd.DataFrame, genes: list[str]) -> float:
    """Mean of the upper-triangle Spearman correlation matrix among
    `genes`. NaN if fewer than 2 genes."""
    if len(genes) < 2:
        return float("nan")
    corr = expr_df[genes].corr(method="spearman").to_numpy()
    upper = corr[np.triu_indices(len(genes), k=1)]
    return float(np.nanmean(upper))


def compute_module_coherence(
    expr_df: pd.DataFrame,
    genes: list[str],
    gene_pool: list[str],
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
) -> dict:
    """Given a seeded `rng`: observed module coherence plus an
    empirical p-value against a random-equal-size-gene-set null drawn
    from `gene_pool`."""
    observed = compute_mean_pairwise_correlation(expr_df, genes)
    pool = np.array([g for g in gene_pool if g not in genes])
    null_values = np.empty(n_permutations)
    for i in range(n_permutations):
        random_genes = rng.choice(pool, size=len(genes), replace=False).tolist()
        null_values[i] = compute_mean_pairwise_correlation(expr_df, random_genes)
    pvalue = float((np.sum(null_values >= observed) + 1) / (n_permutations + 1))
    return {
        "observed_coherence": observed,
        "null_mean": float(np.nanmean(null_values)),
        "pvalue": pvalue,
    }


def build_program_transfer_test(project_root: Path) -> dict:
    gse103322_path = (
        project_root / "data" / "external" / "GSE103322" / "GSE103322_HNSCC_all_data.txt.gz"
    )
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    t_cell_states_path = project_root / "data" / "derived" / "t_cell_states.parquet"
    output_path = project_root / "data" / "derived" / "program_transfer_results.parquet"

    for path, phase in [
        (gse103322_path, None),
        (
            matrix_path,
            "`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`",
        ),
        (t_cell_states_path, "`06_cell_type_annotation/04_resolve_t_cell_substates.R`"),
    ]:
        if not path.exists():
            raise PipelineError(
                f"'{path}' not found."
                + (
                    f" Run {phase} first."
                    if phase
                    else " Run `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py` first."
                )
            )

    # A sizeable (400-gene) random background pool, sampled from
    # GSE103322's gene list, is fetched in addition to the required
    # marker genes -- the null model needs a representative background,
    # not just the handful of marker genes themselves (an earlier
    # version of this module used only the ~29 required genes as its
    # "pool," which would draw null gene sets that overlap heavily with
    # the very programs being tested -- caught and fixed before running).
    all_reference_genes = list_gse103322_genes(gse103322_path)
    background_rng = np.random.default_rng(RNG_SEED)
    candidates = [g for g in all_reference_genes if g not in required_gse103322_genes()]
    background_genes = set(
        background_rng.choice(candidates, size=N_BACKGROUND_GENES, replace=False).tolist()
    )

    reference = parse_gse103322_t_cells(gse103322_path, extra_genes=background_genes)

    # `t_cell_states.parquet`'s index is a plain positional RangeIndex
    # (not cell IDs) -- the cell ID is the `cell_id` column. Using
    # `.index` here previously matched zero cells against
    # `adata.obs_names`, and the first run of this script "succeeded"
    # while actually computing the entire project-side coherence test on
    # an empty dataframe -- a correctness bug, not a cosmetic one; only
    # surfaced because a later fix (`filter_expressed_genes`) turned the
    # resulting empty gene pool into a hard crash instead of a silent
    # wrong answer.
    t_cell_ids = pd.read_parquet(t_cell_states_path)["cell_id"]
    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    gene_pool_all = adata.var_names[adata.var["is_exposure_gene"]].tolist()

    project_t_cells = adata[adata.obs_names.isin(t_cell_ids), gene_pool_all]
    project_expr = pd.DataFrame(
        (
            project_t_cells.layers[layer].toarray()
            if hasattr(project_t_cells.layers[layer], "toarray")
            else project_t_cells.layers[layer]
        ),
        index=project_t_cells.obs_names,
        columns=gene_pool_all,
    )

    reference_gene_pool_raw = [
        g
        for g in reference.columns
        if g
        not in {
            "processed by Maxima enzyme",
            "Lymph node",
            "classified  as cancer cell",
            "classified as non-cancer cells",
            "non-cancer cell type",
        }
    ]
    reference_gene_pool = filter_expressed_genes(reference, reference_gene_pool_raw)
    project_gene_pool = filter_expressed_genes(project_expr, gene_pool_all)

    rows = []
    for program in PROGRAMS_TESTED:
        genes = [g for g in PROGRAM_GENE_SETS[program] if g in gene_pool_all]
        project_result = compute_module_coherence(
            project_expr, genes, project_gene_pool, np.random.default_rng(RNG_SEED)
        )
        reference_genes = [g for g in genes if g in reference.columns]
        reference_result = compute_module_coherence(
            reference, reference_genes, reference_gene_pool, np.random.default_rng(RNG_SEED)
        )

        rows.append(
            {
                "program": program,
                "n_genes": len(genes),
                "project_observed_coherence": project_result["observed_coherence"],
                "project_null_mean": project_result["null_mean"],
                "project_pvalue": project_result["pvalue"],
                "reference_observed_coherence": reference_result["observed_coherence"],
                "reference_null_mean": reference_result["null_mean"],
                "reference_pvalue": reference_result["pvalue"],
                "transfers": bool(
                    project_result["pvalue"] < 0.05 and reference_result["pvalue"] < 0.05
                ),
            }
        )

    # Treg is tested separately (its markers are not in PROGRAM_GENE_SETS).
    treg_genes = [g for g in TREG_MARKERS if g in gene_pool_all]
    project_treg = compute_module_coherence(
        project_expr, treg_genes, project_gene_pool, np.random.default_rng(RNG_SEED)
    )
    reference_treg_genes = [g for g in treg_genes if g in reference.columns]
    reference_treg = compute_module_coherence(
        reference, reference_treg_genes, reference_gene_pool, np.random.default_rng(RNG_SEED)
    )
    rows.append(
        {
            "program": "treg",
            "n_genes": len(treg_genes),
            "project_observed_coherence": project_treg["observed_coherence"],
            "project_null_mean": project_treg["null_mean"],
            "project_pvalue": project_treg["pvalue"],
            "reference_observed_coherence": reference_treg["observed_coherence"],
            "reference_null_mean": reference_treg["null_mean"],
            "reference_pvalue": reference_treg["pvalue"],
            "transfers": bool(project_treg["pvalue"] < 0.05 and reference_treg["pvalue"] < 0.05),
        }
    )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_programs_tested": len(result),
        "n_programs_transferring": int(result["transfers"].sum()),
        "programs_not_transferring": result.loc[~result["transfers"], "program"].tolist(),
        "output_path": str(output_path),
    }
