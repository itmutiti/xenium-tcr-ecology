"""Spatially constrained ligand-receptor interaction scoring (Phase
14.02).

For each panel-complete LR pair (`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`) and each
sender-receiver lineage pair (`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`) whose `relevant_programs`
overlaps that LR pair's `programs` -- a biologically motivated
restriction of the combinatorial space (4 sender-receiver pairs x 9
complete LR pairs = 36 possible combinations naively; program-overlap
restricts this to ~14, avoiding an uninterpretable exhaustive-search
design) -- scores interaction support only across graph-connected
sender-receiver cell pairs (`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated primary graph), not all
possible cell pairs in a section.

**Degree-preserving null, reusing `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s already-calibrated
degree-stratification methodology, not reinventing it:**
`xenium_tcr_ecology.graphs.null_models.compute_degree_strata` is called
directly to bin every cell into local-density quartiles (the same
30um-graph-degree stratification already validated in Phase 9.08's
Type I error / power calibration). The permutation test statistic here
is new (mean ligand(sender) x receptor(receiver) expression product
across sender->receiver graph edges, not `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s distance-
to-nearest-tumour-cell statistic), because this milestone asks a
different question -- but the validated degree-matching concept that
keeps the null valid (preserving each mask's local-density profile, not
just its raw count) is reused exactly as established.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from xenium_tcr_ecology.graphs.null_models import compute_degree_strata
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed

N_PERMUTATIONS = 199
RNG_SEED = get_default_seed()
DEGREE_STRATA_RADIUS_UM = 30.0


def compute_edge_interaction_scores(
    graph: sparse.csr_matrix,
    sender_mask: np.ndarray,
    receiver_mask: np.ndarray,
    ligand_expr: np.ndarray,
    receptor_expr: np.ndarray,
) -> np.ndarray:
    """Pure, testable: per-edge ligand(sender) x receptor(receiver)
    score for every graph edge connecting a sender cell to a receiver
    cell. The graph's symmetric storage (both (i,j) and (j,i) explicit
    for every undirected edge, confirmed in Phase 9.03) means checking
    `row=sender, col=receiver` across all entries already captures both
    physical orientations of each edge."""
    coo = graph.tocoo()
    is_valid_edge = sender_mask[coo.row] & receiver_mask[coo.col]
    return ligand_expr[coo.row[is_valid_edge]] * receptor_expr[coo.col[is_valid_edge]]


def compute_degree_preserving_lr_pvalue(
    graph: sparse.csr_matrix,
    sender_mask: np.ndarray,
    receiver_mask: np.ndarray,
    ligand_expr: np.ndarray,
    receptor_expr: np.ndarray,
    degree_strata: np.ndarray,
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
) -> dict:
    """Pure, testable (given a seeded `rng`): observed mean edge score,
    degree-preserving null distribution (permuting sender and receiver
    labels independently within degree strata), and a one-sided
    empirical p-value (observed unusually high)."""
    n = graph.shape[0]
    observed_scores = compute_edge_interaction_scores(
        graph, sender_mask, receiver_mask, ligand_expr, receptor_expr
    )
    observed = float(observed_scores.mean()) if len(observed_scores) > 0 else float("nan")

    sender_strata_counts = pd.Series(degree_strata[sender_mask]).value_counts()
    receiver_strata_counts = pd.Series(degree_strata[receiver_mask]).value_counts()

    def _permute_mask(strata_counts: pd.Series) -> np.ndarray:
        mask = np.zeros(n, dtype=bool)
        for stratum, n_needed in strata_counts.items():
            pool = np.where(degree_strata == stratum)[0]
            chosen = rng.choice(pool, size=min(int(n_needed), len(pool)), replace=False)
            mask[chosen] = True
        return mask

    null_stats = np.empty(n_permutations)
    for i in range(n_permutations):
        perm_sender = _permute_mask(sender_strata_counts)
        perm_receiver = _permute_mask(receiver_strata_counts)
        perm_scores = compute_edge_interaction_scores(
            graph, perm_sender, perm_receiver, ligand_expr, receptor_expr
        )
        null_stats[i] = float(perm_scores.mean()) if len(perm_scores) > 0 else 0.0

    if np.isnan(observed):
        pvalue = float("nan")
    else:
        pvalue = float((np.sum(null_stats >= observed) + 1) / (n_permutations + 1))

    return {
        "n_edges": int(len(observed_scores)),
        "observed_mean_score": observed,
        "null_mean": float(np.nanmean(null_stats)),
        "pvalue": pvalue,
    }


def find_program_overlap_combinations(
    sender_receiver_pairs: list[dict], lr_pairs: pd.DataFrame
) -> list[dict]:
    """Pure, testable: biologically-motivated combinations of
    (sender-receiver pair, LR pair) where the two `programs` sets
    overlap -- avoids an uninterpretable exhaustive search."""
    combinations = []
    for sr_pair in sender_receiver_pairs:
        sr_programs = set(sr_pair["relevant_programs"])
        for _, lr_pair in lr_pairs.iterrows():
            if not lr_pair["pair_complete"]:
                continue
            lr_programs = set(str(lr_pair["programs"]).split(";")) if lr_pair["programs"] else set()
            if sr_programs & lr_programs:
                combinations.append(
                    {
                        "sender_receiver_pair_id": sr_pair["pair_id"],
                        "sender": sr_pair["sender"],
                        "receiver": sr_pair["receiver"],
                        "lr_pair_id": lr_pair["pair_id"],
                        "ligand": lr_pair["ligand"],
                        "receptor": lr_pair["receptor"],
                    }
                )
    return combinations


def build_spatial_interaction_scores(project_root: Path) -> dict:
    sender_receiver_config_path = project_root / "config" / "sender_receiver_pairs.yaml"
    lr_pairs_path = project_root / "references" / "panel_supported_interactions.tsv"
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    output_path = project_root / "data" / "derived" / "spatial_interaction_scores.parquet"

    for path, phase in [
        (
            sender_receiver_config_path,
            "`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`",
        ),
        (
            lr_pairs_path,
            "`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`",
        ),
        (
            matrix_path,
            "`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`",
        ),
        (final_annotations_path, "`06_cell_type_annotation/06_integrate_annotation_evidence.py`"),
        (
            primary_graphs_dir,
            "`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    sender_receiver_config = yaml.safe_load(sender_receiver_config_path.read_text())
    lr_pairs = pd.read_csv(lr_pairs_path, sep="\t")
    combinations = find_program_overlap_combinations(
        sender_receiver_config["sender_receiver_pairs"], lr_pairs
    )
    if not combinations:
        raise PipelineError(
            "No program-overlap combinations found between sender-receiver pairs and panel-complete LR pairs."
        )

    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    final_annotations = pd.read_parquet(final_annotations_path)

    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for section_dir in sorted(d for d in primary_graphs_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        node_metadata = pd.read_csv(section_dir / "node_metadata.tsv", sep="\t")
        graph = sparse.load_npz(section_dir / "primary_graph.npz")
        graph_binary = graph.copy()
        graph_binary.data = np.ones_like(graph_binary.data)

        section_lineage = final_annotations.reindex(node_metadata["cell_id"])[
            "final_lineage"
        ].to_numpy()
        section_adata = adata[adata.obs_names.isin(node_metadata["cell_id"])]
        section_adata = section_adata[node_metadata["cell_id"].to_numpy()]
        expr = section_adata.layers[layer]
        expr = expr.toarray() if hasattr(expr, "toarray") else expr
        gene_index = {g: i for i, g in enumerate(section_adata.var_names)}

        degree_strata = compute_degree_strata(
            node_metadata[["x_centroid", "y_centroid"]].to_numpy(), DEGREE_STRATA_RADIUS_UM
        )

        for combo in combinations:
            if combo["ligand"] not in gene_index or combo["receptor"] not in gene_index:
                continue
            sender_mask = section_lineage == combo["sender"]
            receiver_mask = section_lineage == combo["receiver"]
            if sender_mask.sum() < 10 or receiver_mask.sum() < 10:
                continue

            ligand_expr = expr[:, gene_index[combo["ligand"]]]
            receptor_expr = expr[:, gene_index[combo["receptor"]]]

            result = compute_degree_preserving_lr_pvalue(
                graph_binary,
                sender_mask,
                receiver_mask,
                ligand_expr,
                receptor_expr,
                degree_strata,
                rng,
            )
            rows.append(
                {
                    "section_id": section_id,
                    "sender_receiver_pair_id": combo["sender_receiver_pair_id"],
                    "lr_pair_id": combo["lr_pair_id"],
                    "n_sender_cells": int(sender_mask.sum()),
                    "n_receiver_cells": int(receiver_mask.sum()),
                    **result,
                }
            )

    result_df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_parquet(output_path)

    return {
        "n_combinations_tested": len(combinations),
        "n_rows": len(result_df),
        "n_sections": int(result_df["section_id"].nunique()) if len(result_df) else 0,
        "n_significant": int((result_df["pvalue"] < 0.05).sum()) if len(result_df) else 0,
        "output_path": str(output_path),
    }
