"""Tumour region mask construction (`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`).

Converts `07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`'s continuous, within-patient-standardised
`malignancy_score` into spatially coherent hard tumour-region calls, per
the specification's "forms spatially coherent tumour regions from
malignant-cell probabilities and removes isolated false positives."

**Threshold, a judgment call made and documented directly:**
`malignancy_score > 0`, not an arbitrary percentile cut.
`malignancy_score` is a mean of four within-patient z-scored components
(`07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py`), so 0 is a meaningful reference point -- "no net
evidence of malignancy relative to this patient's own non-epithelial
reference baseline" -- not an arbitrary round number. A percentile-based
cut (e.g. `malignancy_probability > 0.5`) was considered and rejected:
`malignancy_probability` is a percentile rank, uniform on [0, 1] by
construction regardless of the underlying score's shape, so any fixed
percentile cut selects that exact fraction of cells always, whether the
true malignant fraction in a given section is 10% or 90% -- it carries no
information about whether a cell's evidence is actually net positive or
net negative, only its rank among other epithelial cells. Checked
against the distribution before finalising: median
`malignancy_score` is -0.019, i.e. close to but not exactly 0, so this is
not merely a relabelled median split.

**Spatial coherence, two stages:**
1. Neighbour-majority smoothing: a cell's smoothed call is the majority
   vote among its own k=10 spatial nearest neighbours' raw hard calls
   (same k and rationale as `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s `spatial_consistency` --
   tissue is spatially organised, so an isolated single-cell disagreement
   with its whole local neighbourhood is more likely segmentation/
   technical noise than genuine fine-grained heterogeneity at this
   scale).
2. Connected-component filtering: smoothed-malignant cells within
   `ADJACENCY_RADIUS_UM` of each other (`04_quality_control/04_estimate_transcript_spillover.py`'s data-grounded
   search radius, ~one median cell diameter) are grouped into discrete
   spatial regions (`scipy.sparse.csgraph.connected_components`).
   Regions smaller than `MIN_REGION_SIZE_CELLS` are excluded as isolated
   false positives -- the specification's explicit requirement -- their
   cells revert to "not part of a tumour region," not deleted from the
   dataset.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors

from xenium_tcr_ecology.infra.exceptions import PipelineError

SPATIAL_K_NEIGHBORS = 10
ADJACENCY_RADIUS_UM = 10.0
MIN_REGION_SIZE_CELLS = 5


def compute_smoothed_calls(
    x: np.ndarray, y: np.ndarray, raw_calls: np.ndarray, k: int = SPATIAL_K_NEIGHBORS
) -> np.ndarray:
    """Majority-vote smoothing of a boolean per-cell call over its own k
    spatial nearest neighbours (excluding itself)."""
    n = len(x)
    if n == 0:
        return np.array([], dtype=bool)
    coords = np.column_stack([x, y])
    k_eff = min(k + 1, n)
    if k_eff < 2:
        return raw_calls.copy()
    nn = NearestNeighbors(n_neighbors=k_eff).fit(coords)
    _, indices = nn.kneighbors(coords)
    neighbor_calls = raw_calls[indices[:, 1:]]
    return neighbor_calls.mean(axis=1) >= 0.5


def label_connected_regions(
    x: np.ndarray, y: np.ndarray, is_malignant: np.ndarray, radius_um: float = ADJACENCY_RADIUS_UM
) -> np.ndarray:
    """Connected-component region ID for cells with `is_malignant == True`
    that are within `radius_um` of at least one other such cell (via a
    radius graph, not a fixed-k graph -- region shape should not be forced
    to a constant node degree). Returns an int array the same length as
    `x`/`y`/`is_malignant`; -1 for cells not part of any malignant
    region."""
    n = len(x)
    region_id = np.full(n, -1, dtype=int)
    malignant_idx = np.where(is_malignant)[0]
    if len(malignant_idx) == 0:
        return region_id

    coords = np.column_stack([x, y])[malignant_idx]
    nn = NearestNeighbors(radius=radius_um).fit(coords)
    adjacency = nn.radius_neighbors_graph(coords, mode="connectivity")
    n_components, labels = connected_components(coo_matrix(adjacency), directed=False)
    region_id[malignant_idx] = labels
    return region_id


def filter_small_regions(
    region_id: np.ndarray, min_size: int = MIN_REGION_SIZE_CELLS
) -> np.ndarray:
    """Reassigns cells in regions smaller than `min_size` back to -1
    ("not part of a tumour region") -- the specification's explicit
    "removes isolated false positives" step. Region IDs are not
    renumbered to stay contiguous after filtering; only membership is
    affected."""
    result = region_id.copy()
    valid = result >= 0
    if not valid.any():
        return result
    sizes = pd.Series(result[valid]).value_counts()
    small_regions = sizes.index[sizes < min_size]
    result[np.isin(result, small_regions)] = -1
    return result


def build_tumour_region_masks(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    malignancy_scores_path = project_root / "data" / "derived" / "malignancy_scores.parquet"
    output_dir = project_root / "data" / "derived" / "tumour_masks"

    for p in (matrix_path, malignancy_scores_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    scores = pd.read_parquet(malignancy_scores_path)
    scores = scores.join(adata.obs[["x_centroid", "y_centroid"]])

    output_dir.mkdir(parents=True, exist_ok=True)
    section_rows = []
    section_outputs = []
    for section_id, group in scores.groupby("section_id", observed=True):
        x = group["x_centroid"].to_numpy()
        y = group["y_centroid"].to_numpy()
        raw_call = (group["malignancy_score"] > 0).to_numpy()

        smoothed_call = compute_smoothed_calls(x, y, raw_call)
        region_id = label_connected_regions(x, y, smoothed_call)
        region_id = filter_small_regions(region_id)

        section_result = pd.DataFrame(
            {
                "raw_malignant_call": raw_call,
                "smoothed_malignant_call": smoothed_call,
                "region_id": region_id,
                "in_tumour_region": region_id >= 0,
            },
            index=group.index,
        )
        section_result.to_parquet(output_dir / f"{section_id}.parquet")
        section_outputs.append(section_result)

        n_regions = len(np.unique(region_id[region_id >= 0]))
        section_rows.append(
            {
                "section_id": section_id,
                "n_cells": len(group),
                "n_raw_malignant": int(raw_call.sum()),
                "n_smoothed_malignant": int(smoothed_call.sum()),
                "n_in_tumour_region": int((region_id >= 0).sum()),
                "n_regions": n_regions,
                "n_isolated_removed": int(smoothed_call.sum() - (region_id >= 0).sum()),
            }
        )

    summary_df = pd.DataFrame(section_rows)
    summary_path = output_dir / "_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    return {
        "n_sections": len(summary_df),
        "n_cells_total": int(summary_df["n_cells"].sum()),
        "n_in_tumour_region_total": int(summary_df["n_in_tumour_region"].sum()),
        "fraction_in_tumour_region": round(
            float(summary_df["n_in_tumour_region"].sum() / summary_df["n_cells"].sum()), 4
        ),
        "n_regions_total": int(summary_df["n_regions"].sum()),
        "n_isolated_removed_total": int(summary_df["n_isolated_removed"].sum()),
        "output_dir": str(output_dir),
        "summary_path": str(summary_path),
    }
