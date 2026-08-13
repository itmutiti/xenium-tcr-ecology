"""Programme activity specifically at fibroblast/suppressive-myeloid
barrier interfaces, contrasting "excluded" vs "engaged" T-cell clones
(`14_spatial_interactions_and_barriers/04_analyse_barrier_pathways.py`).

Exploratory -- not a prespecified confirmatory analysis (only Phase
14.03's `q3_barrier_topology_confirmatory` is prespecified; see
governance/analysis_registry.tsv). Follows from `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s
significant finding that `suppressive_myeloid_barrier_fraction`
predicts lower clone-tumour engagement: this milestone asks whether
the same barrier cells (fibroblast/suppressive-myeloid lineage cells
lying on a clone's shortest graph path to the nearest tumour cell,
`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s method) differ in checkpoint/
chemokine/interferon/antigen-presentation programme activity depending
on whether they interpose an excluded clone (engagement_ratio below
the cohort median) or an engaged clone (above median) from the tumour.

**Programme set, matching `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s scaffold wording
("checkpoint, chemokine, TGF-beta, interferon and antigen-presentation
programs"):** `exhaustion_score` (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s name for the
checkpoint programme), `interferon_score`, `antigen_presentation_score`
are reused directly from `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s already-computed `data/derived/
program_scores.parquet` (confirmed 1:1 cell_id coverage against
`primary_analysis_matrix.h5ad` before joining -- no recomputation).
`chemokine_score` did not previously exist as a per-cell score (Phase
14.00 only defined `CHEMOKINE_GENE_SET` as a gene list) -- computed here
for the first time, using the same `scanpy.tl.score_genes` /
`is_exposure_gene`-pool technique `05_preprocessing_and_normalisation/03_calculate_program_scores.py` itself used (not the
GSE103322 reference's simplified z-score proxy,
`external_checkpoint.bulk_reference.compute_module_score` -- that simplification
was scoped to the bulk-reference cross-platform comparison and would
not be a like-for-like method for this project's own data, which
already has the full panel available for a control-matched score).
`tgf_beta` is skipped explicitly (`TGFB_GENE_SET` is empty,
`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`), not silently omitted.

**Interface-cell identification, reusing `11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s already-
validated method rather than reimplementing it:** `assign_barrier_group`,
`trace_intermediate_path`, and the lineage/substate constants are
imported directly from `barrier_metrics.py`. The per-section
multi-source Dijkstra solve is repeated here (`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`'s
`build_clone_barrier_metrics` did not persist predecessor arrays or
interface cell identities, only aggregate barrier fractions) -- a
necessary recomputation of the same, already-validated method, not a
new or divergent one.

**Pseudoreplication decision:** a single barrier cell can lie on the
shortest path for multiple clones (a "hub" barrier cell). Rows are
deduplicated to one (cell_id, section_id, clone_class) tuple before the
group comparison -- a cell bordering several excluded clones is counted
once in the "excluded" group, not once per referencing clone -- so the
group comparison is not inflated by a small number of highly-connected
barrier cells.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
from scipy.stats import mannwhitneyu

from xenium_tcr_ecology.clone_ecology.barrier_metrics import (
    TUMOUR_LINEAGE,
    assign_barrier_group,
    trace_intermediate_path,
)
from xenium_tcr_ecology.clone_ecology.spatial_descriptors import filter_to_primary_cohort
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed
from xenium_tcr_ecology.interactions.sender_receiver_pairs import CHEMOKINE_GENE_SET, TGFB_GENE_SET

TARGET_BARRIER_GROUPS = ("fibroblast", "suppressive_myeloid")
PROGRAM_SCORE_COLUMNS = {
    "checkpoint": "exhaustion_score",
    "interferon": "interferon_score",
    "antigen_presentation": "antigen_presentation_score",
    "chemokine": "chemokine_score",
}
RNG_SEED = get_default_seed()


def classify_clone_engagement(engagement_ratio: pd.Series) -> pd.Series:
    """Median-split classification -- values below the cohort median
    are "excluded", values at or above the median are "engaged" (an
    explicit boundary convention, not left ambiguous at the median
    value itself)."""
    median = engagement_ratio.median()
    return pd.Series(
        np.where(engagement_ratio < median, "excluded", "engaged"), index=engagement_ratio.index
    )


def identify_barrier_interface_cells(
    predecessors: np.ndarray,
    tumour_idx_set: set[int],
    clone_cell_indices: list[int],
    barrier_group: np.ndarray,
    target_groups: tuple[str, ...] = TARGET_BARRIER_GROUPS,
) -> set[int]:
    """Union of interface-cell indices across every clone cell's
    shortest path to the tumour, restricted to `target_groups`
    barrier-lineage cells only."""
    interface_cells: set[int] = set()
    for cell_idx in clone_cell_indices:
        path_idx = trace_intermediate_path(predecessors, cell_idx, tumour_idx_set)
        for idx in path_idx:
            if barrier_group[idx] in target_groups:
                interface_cells.add(idx)
    return interface_cells


def _benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Pure, testable: standard BH step-up FDR adjustment."""
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order] * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty(n)
    result[order] = adjusted
    return result


def compare_interface_programs(
    interface_table: pd.DataFrame, program_scores: pd.DataFrame, programs: list[str]
) -> pd.DataFrame:
    """Two-sided Mann-Whitney U comparison of excluded- vs.
    engaged-clone interface cells for each programme score column,
    with BH FDR adjustment across `programs`."""
    deduped = interface_table.drop_duplicates(subset=["cell_id", "clone_class"])
    merged = deduped.merge(program_scores, left_on="cell_id", right_index=True, how="inner")

    rows = []
    for program in programs:
        excluded_scores = merged.loc[merged["clone_class"] == "excluded", program].dropna()
        engaged_scores = merged.loc[merged["clone_class"] == "engaged", program].dropna()
        if len(excluded_scores) >= 2 and len(engaged_scores) >= 2:
            _, pvalue = mannwhitneyu(excluded_scores, engaged_scores, alternative="two-sided")
        else:
            pvalue = float("nan")
        rows.append(
            {
                "program": program,
                "n_excluded_interface_cells": len(excluded_scores),
                "n_engaged_interface_cells": len(engaged_scores),
                "mean_excluded": (
                    float(excluded_scores.mean()) if len(excluded_scores) else float("nan")
                ),
                "mean_engaged": (
                    float(engaged_scores.mean()) if len(engaged_scores) else float("nan")
                ),
                "pvalue": pvalue,
            }
        )
    result = pd.DataFrame(rows)
    result["pvalue_bh"] = _benjamini_hochberg(result["pvalue"].to_numpy())
    return result


def build_barrier_pathways(project_root: Path) -> dict:
    engagement_path = project_root / "data" / "derived" / "clone_tumour_engagement.parquet"
    barrier_path = project_root / "data" / "derived" / "clone_barrier_metrics.parquet"
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    program_scores_path = project_root / "data" / "derived" / "program_scores.parquet"
    output_path = project_root / "data" / "derived" / "barrier_pathways.parquet"

    for path, phase in [
        (engagement_path, "`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`"),
        (
            barrier_path,
            "`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`",
        ),
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
        (sample_manifest_path, None),
        (final_annotations_path, "`06_cell_type_annotation/06_integrate_annotation_evidence.py`"),
        (
            primary_graphs_dir,
            "`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`",
        ),
        (
            matrix_path,
            "`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`",
        ),
        (
            program_scores_path,
            "`05_preprocessing_and_normalisation/03_calculate_program_scores.py`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found." + (f" Run {phase} first." if phase else ""))

    if len(TGFB_GENE_SET) == 0:
        tgf_beta_note = "skipped: TGFB_GENE_SET has zero panel coverage (`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`)"
    else:
        tgf_beta_note = None

    engagement = pd.read_parquet(engagement_path)
    barrier = pd.read_parquet(barrier_path)
    clone_data = engagement.merge(
        barrier[
            [
                "clone_id",
                "section_id",
                "fibroblast_barrier_fraction",
                "suppressive_myeloid_barrier_fraction",
            ]
        ],
        on=["clone_id", "section_id"],
    )
    clone_data = clone_data.dropna(
        subset=["fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction"]
    )
    if len(clone_data) == 0:
        raise PipelineError(
            "No clone-section rows with a defined barrier fraction -- cannot classify interfaces."
        )
    clone_data["clone_class"] = classify_clone_engagement(clone_data["engagement_ratio"])
    clone_class_lookup = clone_data.set_index(["clone_id", "section_id"])["clone_class"]

    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
    primary_sections = set(filter_to_primary_cohort(sample_manifest)["section_id"])

    resolved = pd.read_parquet(resolved_calls_path)
    clonal = resolved[resolved["resolution"].isin(["singlet", "low_confidence"])].copy()
    clonal["clone_id"] = clonal["detected_probes"]
    high_confidence_clones = pd.read_parquet(high_confidence_clones_path)
    clonal = clonal[
        clonal["clone_id"].isin(set(high_confidence_clones["clone_id"]))
        & clonal["section_id"].isin(primary_sections)
    ]

    final_annotations = pd.read_parquet(final_annotations_path)

    adata = ad.read_h5ad(matrix_path)
    layer = adata.uns["primary_normalization_layer"]
    gene_pool = adata.var_names[adata.var["is_exposure_gene"]].tolist()
    sc.tl.score_genes(
        adata,
        gene_list=CHEMOKINE_GENE_SET,
        gene_pool=gene_pool,
        layer=layer,
        score_name="chemokine_score",
        random_state=RNG_SEED,
    )
    chemokine_scores = adata.obs["chemokine_score"]

    existing_program_scores = pd.read_parquet(program_scores_path)
    program_scores = existing_program_scores.join(chemokine_scores)

    interface_rows = []
    for section_id in sorted(clonal["section_id"].unique()):
        node_metadata = pd.read_csv(primary_graphs_dir / section_id / "node_metadata.tsv", sep="\t")
        graph = sparse.load_npz(primary_graphs_dir / section_id / "primary_graph.npz")

        section_annotations = final_annotations.reindex(node_metadata["cell_id"])
        lineage = section_annotations["final_lineage"].to_numpy()
        substate = section_annotations["final_substate"].to_numpy()
        barrier_group = assign_barrier_group(lineage, substate)

        tumour_idx = np.flatnonzero(lineage == TUMOUR_LINEAGE)
        if len(tumour_idx) == 0:
            continue
        tumour_idx_set = set(tumour_idx.tolist())

        _, predecessors, _ = dijkstra(
            graph, directed=False, indices=tumour_idx, min_only=True, return_predecessors=True
        )

        cell_id_to_idx = {cid: i for i, cid in enumerate(node_metadata["cell_id"])}
        section_clonal = clonal[clonal["section_id"] == section_id]

        for clone_id, group in section_clonal.groupby("clone_id", observed=True):
            clone_class = clone_class_lookup.get((clone_id, section_id))
            if clone_class is None:
                continue
            cell_indices = [cell_id_to_idx[cid] for cid in group.index if cid in cell_id_to_idx]
            if len(cell_indices) == 0:
                continue
            interface_indices = identify_barrier_interface_cells(
                predecessors, tumour_idx_set, cell_indices, barrier_group
            )
            for idx in interface_indices:
                interface_rows.append(
                    {
                        "cell_id": node_metadata["cell_id"].iloc[idx],
                        "section_id": section_id,
                        "clone_id": clone_id,
                        "clone_class": clone_class,
                    }
                )

    interface_table = pd.DataFrame(interface_rows)
    if len(interface_table) == 0:
        raise PipelineError(
            "No fibroblast/suppressive-myeloid interface cells found for any classified clone."
        )

    programs_to_test = [c for c in PROGRAM_SCORE_COLUMNS.values()]
    result = compare_interface_programs(interface_table, program_scores, programs_to_test)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_classified_clone_sections": len(clone_data),
        "n_interface_cell_clone_rows": len(interface_table),
        "n_unique_interface_cells": interface_table["cell_id"].nunique(),
        "n_programs_tested": len(result),
        "n_significant_bh": int((result["pvalue_bh"] < 0.05).sum()),
        "tgf_beta_note": tgf_beta_note,
        "output_path": str(output_path),
    }
