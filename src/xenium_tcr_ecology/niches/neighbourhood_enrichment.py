"""Cell-type neighbourhood enrichment (`10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py`).

Tests pairwise cell-type adjacency (`06_cell_type_annotation/06_integrate_annotation_evidence.py`'s `final_lineage`, all 12
major lineages) using within-section constrained permutations -- the
exact null model methodology calibrated in `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py` (`constrained
permutation`: labels randomly reassigned among graph nodes, unconstrained
beyond class counts). Rather than reimplementing this by hand for all 66
lineage pairs across 18 sections, this module uses `squidpy.gr.
nhood_enrichment` directly on `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s calibrated primary graph
(30um, gap-pruned, patient-separated) -- the same established library
already used and validated for spatial permutation testing
elsewhere in this project (`07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`, 8.04), and methodologically
equivalent to `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s own "constrained permutation" null (unrestricted
label reassignment among graph nodes, same graph, same style of
permutation test).

**Within-section, not pooled:** run independently per section (never
pooling graphs across sections, since `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s primary graph already
never connects cells across sections) -- satisfies the blueprint's own
"within-section" qualifier structurally, not via a separate constraint
parameter.

**Aggregation across sections:** for each of the 66 unique lineage pairs,
the per-section z-scores are summarised (median, and the fraction of
sections agreeing in sign) across all sections where both lineages have
at least `MIN_CELLS_PER_LINEAGE` cells -- a pair's enrichment/depletion
call is only as trustworthy as its consistency across the separate
tissue sections it is measured in, not a single pooled number.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

N_PERMUTATIONS = 199
MIN_CELLS_PER_LINEAGE = 10


def build_section_adata(node_metadata: pd.DataFrame, graph: sparse.csr_matrix):
    import anndata as ad

    adata = ad.AnnData(
        X=np.zeros((len(node_metadata), 1), dtype=np.float32), obs=node_metadata.copy()
    )
    adata.obs["final_lineage"] = adata.obs["final_lineage"].astype("category")
    adata.obsp["spatial_connectivities"] = (graph > 0).astype(int)
    return adata


def compute_section_enrichment(
    node_metadata: pd.DataFrame,
    graph: sparse.csr_matrix,
    n_perms: int = N_PERMUTATIONS,
    seed: int = 0,
) -> pd.DataFrame:
    """Pure(ish), testable per-section enrichment computation. Returns a
    long-format DataFrame (lineage_a, lineage_b, zscore, count) for
    every pair with both lineages meeting MIN_CELLS_PER_LINEAGE."""
    import squidpy as sq

    lineage_counts = node_metadata["final_lineage"].value_counts()
    eligible_lineages = lineage_counts[lineage_counts >= MIN_CELLS_PER_LINEAGE].index.tolist()
    if len(eligible_lineages) < 2:
        return pd.DataFrame(columns=["lineage_a", "lineage_b", "zscore", "count"])

    mask = node_metadata["final_lineage"].isin(eligible_lineages).to_numpy()
    sub_metadata = node_metadata.loc[mask]
    sub_graph = graph[mask][:, mask]

    adata = build_section_adata(sub_metadata, sub_graph)
    result = sq.gr.nhood_enrichment(
        adata,
        cluster_key="final_lineage",
        n_perms=n_perms,
        seed=seed,
        copy=True,
        show_progress_bar=False,
    )
    categories = adata.obs["final_lineage"].cat.categories.tolist()

    rows = []
    for i, lineage_a in enumerate(categories):
        for j, lineage_b in enumerate(categories):
            if j < i:
                continue
            rows.append(
                {
                    "lineage_a": lineage_a,
                    "lineage_b": lineage_b,
                    "zscore": float(result.zscore[i, j]),
                    "count": int(result.counts[i, j]),
                }
            )
    return pd.DataFrame(rows)


def build_neighbourhood_enrichment(project_root: Path) -> dict:
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    output_path = project_root / "data" / "derived" / "neighbourhood_enrichment.parquet"

    if not primary_graphs_dir.is_dir():
        raise PipelineError(
            f"'{primary_graphs_dir}' not found. Run `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py` first."
        )

    all_rows = []
    for section_dir in sorted(d for d in primary_graphs_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        node_metadata = pd.read_csv(section_dir / "node_metadata.tsv", sep="\t", index_col=0)
        graph = sparse.load_npz(section_dir / "primary_graph.npz")

        section_result = compute_section_enrichment(node_metadata, graph)
        section_result["section_id"] = section_id
        all_rows.append(section_result)

    per_section = pd.concat(all_rows, ignore_index=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    per_section.to_parquet(output_path)

    summary = (
        per_section.groupby(["lineage_a", "lineage_b"])["zscore"]
        .agg(median_zscore="median", n_sections="count", frac_positive=lambda s: (s > 0).mean())
        .reset_index()
    )
    summary_path = project_root / "data" / "derived" / "neighbourhood_enrichment_summary.parquet"
    summary.to_parquet(summary_path)

    top_enriched = (
        summary[summary["lineage_a"] != summary["lineage_b"]]
        .sort_values("median_zscore", ascending=False)
        .head(5)
    )
    top_depleted = (
        summary[summary["lineage_a"] != summary["lineage_b"]].sort_values("median_zscore").head(5)
    )

    return {
        "n_sections": per_section["section_id"].nunique(),
        "n_lineage_pairs": len(summary),
        "top_enriched_pairs": top_enriched[["lineage_a", "lineage_b", "median_zscore"]].to_dict(
            "records"
        ),
        "top_depleted_pairs": top_depleted[["lineage_a", "lineage_b", "median_zscore"]].to_dict(
            "records"
        ),
        "output_path": str(output_path),
        "summary_path": str(summary_path),
    }
