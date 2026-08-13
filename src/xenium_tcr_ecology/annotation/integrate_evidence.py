"""Integration of marker, reference, cluster and spatial evidence into
final cell annotations with calibrated confidence and ambiguity (`06_cell_type_annotation/06_integrate_annotation_evidence.py`).

The primary label for every cell is the marker-based major-lineage argmax
(`06_cell_type_annotation/02_score_major_lineages.py`) -- the only evidence source available for all 1,120,780 cells
(reference-transfer is T-cell-only, per `06_cell_type_annotation/03_map_external_scrna_reference.py`; substates are
compartment-gated, per `06_cell_type_annotation/04_resolve_t_cell_substates.R`, `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R`). Three independent, [0,1]-scaled
consistency measures are combined into one confidence score:

  1. marker_margin: normalised gap between the top and second-best lineage
     score (a large margin means the marker call is unambiguous on its own
     terms).
  2. cluster_consensus: fraction of a cell's own Leiden cluster (Phase
     6.01, resolution 0.5, chosen as a moderate, non-extreme resolution)
     that shares its lineage call -- internal consistency with exploratory
     transcriptomic structure.
  3. spatial_consistency: fraction of a cell's k=10 nearest spatial
     neighbours *within the same section* (tissue architecture is not
     comparable across sections) that share its lineage call --
     tissue is spatially organised (tumour nests, immune stroma), so an
     isolated single-cell disagreement with its whole local neighbourhood
     is more likely a technical/segmentation artefact than a genuine
     minority cell type, without treating spatial disagreement as
     automatically wrong (an isolated immune cell inside a tumour nest
     is exactly the kind of biology this project cares about, so
     spatial evidence is one vote among three, not an overriding one).

`confidence = mean(marker_margin, cluster_consensus, spatial_consistency)`.
Cells below `AMBIGUITY_THRESHOLD` are flagged `is_ambiguous`; per the
specification's "fine subtypes are withheld when the panel cannot
support them," ambiguous cells' substate columns are cleared (kept as
NA) even if `06_cell_type_annotation/04_resolve_t_cell_substates.R`, `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R` produced a substate call for them, but the major-lineage
label itself is retained (a cell is still, at minimum, "the lineage its
strongest available evidence points to," not un-labelled).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from xenium_tcr_ecology.infra.exceptions import PipelineError

CLUSTER_RESOLUTION_COL = "joint_leiden_res0.5"
SPATIAL_K_NEIGHBORS = 10
# Recalibrated against the confidence distribution, not chosen a priori:
# an initial value of 0.4 was picked before any data existed to check it
# against, and produced an implausible 57.5% "ambiguous" rate once
# marker_margin's normalisation bug was fixed -- confidence's
# 25th percentile on the data is 0.198, so 0.4 was flagging closer to the
# population's *median* as ambiguous, not a genuine minority. 0.2 flags
# roughly the bottom quartile -- a defensible, documented "uncertain
# minority" cut, consistent with the specification's framing of
# ambiguity as a flag on some cells, not most. cluster_consensus and
# spatial_consistency are not further rescaled to compensate: unlike
# marker_margin's raw score-difference units (which have no natural [0,1]
# interpretation without normalisation), these are already meaningful,
# interpretable absolute fractions (e.g. "70% of my spatial neighbours
# agree with me"), and their skewed-low distribution reflects genuine
# biological/technical structure (Leiden clusters and marker-based
# lineage calls are two different, imperfectly-aligned views of the same
# cells; tissue has genuine spatial heterogeneity at a k=10-neighbour
# scale) -- not a normalisation artefact to paper over.
AMBIGUITY_THRESHOLD = 0.2


def compute_marker_margin(lineage_scores: pd.DataFrame) -> pd.Series:
    lineage_cols = [c for c in lineage_scores.columns if c.endswith("_lineage_score")]
    values = lineage_scores[lineage_cols].to_numpy()
    sorted_vals = np.sort(values, axis=1)
    raw_margin = pd.Series(sorted_vals[:, -1] - sorted_vals[:, -2], index=lineage_scores.index)
    # Percentile-rank, not min-max, scaling to [0, 1]: confirmed on the
    # data that the raw margin distribution is heavily right-skewed (a
    # small number of very-high-margin cells stretch a min-max scale so
    # far that the bulk of the population -- moderate-confidence
    # cells -- gets compressed toward 0, producing a mean marker_margin of
    # 0.071 and driving an implausible 83.2% ambiguity rate on the first
    # run). Percentile rank is robust to this: it depends only on each
    # cell's relative order in the population, not the raw scale of a
    # few extreme outliers, and is naturally uniform-if-monotonic on
    # [0, 1] regardless of the underlying distribution's skew.
    return raw_margin.rank(pct=True)


def compute_cluster_consensus(
    clustering_assignments: pd.DataFrame, argmax_lineage: pd.Series
) -> pd.Series:
    """For each cell, the fraction of its own Leiden cluster that shares
    its own lineage label -- not the cluster's overall majority-agreement
    rate applied uniformly to every member. A cell in the cluster's
    majority label gets a high score; a minority cell in the same cluster
    correctly gets a lower one, since it disagrees with most of its own
    cluster-mates even though the cluster as a whole is mostly consistent.
    """
    if CLUSTER_RESOLUTION_COL not in clustering_assignments.columns:
        raise PipelineError(f"'{CLUSTER_RESOLUTION_COL}' not found in clustering_assignments.")
    df = pd.DataFrame(
        {"cluster": clustering_assignments[CLUSTER_RESOLUTION_COL], "lineage": argmax_lineage}
    )
    pair_counts = df.groupby(["cluster", "lineage"]).size().rename("pair_count")
    cluster_size = df.groupby("cluster").size().rename("cluster_size")
    df = df.join(pair_counts, on=["cluster", "lineage"]).join(cluster_size, on="cluster")
    return df["pair_count"] / df["cluster_size"]


def compute_spatial_consistency(
    section_ids: pd.Series,
    x: np.ndarray,
    y: np.ndarray,
    argmax_lineage: pd.Series,
    k: int = SPATIAL_K_NEIGHBORS,
) -> pd.Series:
    result = pd.Series(index=section_ids.index, dtype=float)
    for section_id in section_ids.unique():
        mask = (section_ids == section_id).to_numpy()
        coords = np.column_stack([x[mask], y[mask]])
        labels = argmax_lineage[mask].to_numpy()
        n = len(coords)
        k_eff = min(k + 1, n)
        if k_eff < 2:
            result.loc[section_ids.index[mask]] = 1.0
            continue
        nn = NearestNeighbors(n_neighbors=k_eff).fit(coords)
        _, indices = nn.kneighbors(coords)
        neighbor_labels = labels[indices[:, 1:]]
        same_label = (neighbor_labels == labels[:, None]).mean(axis=1)
        result.loc[section_ids.index[mask]] = same_label
    return result


def integrate_evidence(
    lineage_scores: pd.DataFrame,
    clustering_assignments: pd.DataFrame,
    section_ids: pd.Series,
    x: np.ndarray,
    y: np.ndarray,
    t_cell_states: pd.DataFrame | None = None,
    tme_substates: pd.DataFrame | None = None,
) -> pd.DataFrame:
    lineage_cols = [c for c in lineage_scores.columns if c.endswith("_lineage_score")]
    if len(lineage_cols) == 0:
        raise PipelineError("lineage_scores has no '*_lineage_score' columns.")

    argmax_lineage = (
        lineage_scores[lineage_cols].idxmax(axis=1).str.replace("_lineage_score", "", regex=False)
    )

    marker_margin = compute_marker_margin(lineage_scores)
    cluster_consensus = compute_cluster_consensus(clustering_assignments, argmax_lineage)
    spatial_consistency = compute_spatial_consistency(section_ids, x, y, argmax_lineage)

    confidence = (marker_margin + cluster_consensus + spatial_consistency) / 3.0
    is_ambiguous = confidence < AMBIGUITY_THRESHOLD

    result = pd.DataFrame(
        {
            "final_lineage": argmax_lineage,
            "marker_margin": marker_margin,
            "cluster_consensus": cluster_consensus,
            "spatial_consistency": spatial_consistency,
            "confidence": confidence,
            "is_ambiguous": is_ambiguous,
        },
        index=lineage_scores.index,
    )

    result["final_substate"] = pd.Series(pd.NA, index=result.index, dtype="object")
    if t_cell_states is not None:
        result.loc[t_cell_states.index, "final_substate"] = t_cell_states["t_cell_state"]
    if tme_substates is not None:
        result.loc[tme_substates.index, "final_substate"] = tme_substates["tme_substate"]
    result.loc[result["is_ambiguous"], "final_substate"] = pd.NA

    return result


def build_final_annotations_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    lineage_scores_path = project_root / "data" / "derived" / "lineage_scores.parquet"
    clustering_path = project_root / "data" / "derived" / "clustering_assignments.parquet"
    t_cell_states_path = project_root / "data" / "derived" / "t_cell_states.parquet"
    tme_substates_path = project_root / "data" / "derived" / "tme_substates.parquet"
    output_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"

    for p in (
        matrix_path,
        lineage_scores_path,
        clustering_path,
        t_cell_states_path,
        tme_substates_path,
    ):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    adata = ad.read_h5ad(matrix_path)
    lineage_scores = pd.read_parquet(lineage_scores_path)
    clustering_assignments = pd.read_parquet(clustering_path)
    t_cell_states = pd.read_parquet(t_cell_states_path).set_index("cell_id")
    tme_substates = pd.read_parquet(tme_substates_path).set_index("cell_id")

    section_ids = adata.obs["section_id"]
    x = adata.obs["x_centroid"].to_numpy()
    y = adata.obs["y_centroid"].to_numpy()

    result = integrate_evidence(
        lineage_scores, clustering_assignments, section_ids, x, y, t_cell_states, tme_substates
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_cells": len(result),
        "n_ambiguous": int(result["is_ambiguous"].sum()),
        "fraction_ambiguous": round(float(result["is_ambiguous"].mean()), 4),
        "final_lineage_counts": result["final_lineage"].value_counts().to_dict(),
        "mean_confidence": round(float(result["confidence"].mean()), 4),
        "n_with_substate": int(result["final_substate"].notna().sum()),
    }
