"""Tumour boundary geometry and signed distance-to-boundary (`07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`).

Builds an actual polygon geometry for each section's tumour mass from
`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`'s per-cell region masks, then computes a signed distance-to-
boundary value for every cell in the section -- all lineages, not only
epithelial cells (the specification's "for all cells" is taken literally:
a T cell's distance to the nearest tumour boundary is exactly the
downstream measurement this whole project's clone-tumour-engagement
question needs, `11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`, Clone Ecology Confirmatory Models).

**Boundary construction:** each in-tumour-region cell (`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`) is
buffered by `BOUNDARY_BUFFER_UM` (half of `04_quality_control/04_estimate_transcript_spillover.py`'s data-
grounded adjacency radius, i.e. ~one median cell radius rather than one
cell diameter -- deliberately smaller than the 10um radius used to define
region connectivity in `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`, so the resulting mask traces the
cell footprints without artificially inflating the tumour mass by the
full inter-cell gap distance) and the buffered disks are unioned
(`shapely.unary_union`) into one (Multi)Polygon per section -- a simple,
robust way to turn a scattered point set into a spatially coherent area,
avoiding the fragility of alpha-shape/concave-hull algorithms for this
purpose.

**Signed distance convention, stated explicitly:** distance to the
tumour-mask boundary line (`polygon.boundary`, not the filled area -- the
filled area gives 0 for every interior point, which is not a distance-to-
boundary at all), signed negative for cells inside the tumour mask and
positive for cells outside -- the standard signed-distance-field
convention. `distance_to_tumour_boundary_um` is always the unsigned
magnitude; `signed_distance_to_tumour_boundary_um` carries the sign;
`is_inside_tumour_region` is the direct boolean for callers who want it
without inferring it from the sign.

Sections with no tumour regions at all (`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`'s result for
P01_run1: zero regions survived size-filtering) get every cell marked
`is_inside_tumour_region = False` and a null distance -- there is no
tumour boundary to measure distance to in that section, which is a
correct, not missing, result and is not imputed or dropped silently.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import shapely

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.region_masks import ADJACENCY_RADIUS_UM

BOUNDARY_BUFFER_UM = ADJACENCY_RADIUS_UM / 2.0


def build_tumour_mask_polygon(
    in_region_x: np.ndarray, in_region_y: np.ndarray, buffer_um: float = BOUNDARY_BUFFER_UM
):
    """Returns a single (Multi)Polygon covering all in-tumour-region cells
    in one section, or None if there are no such cells."""
    if len(in_region_x) == 0:
        return None
    points = shapely.points(in_region_x, in_region_y)
    disks = shapely.buffer(points, buffer_um)
    return shapely.unary_union(disks)


def compute_signed_distances(x: np.ndarray, y: np.ndarray, tumour_mask_polygon) -> pd.DataFrame:
    """Pure, testable signed-distance computation for one section's worth
    of cells against an already-built tumour mask polygon (or None)."""
    n = len(x)
    if tumour_mask_polygon is None or tumour_mask_polygon.is_empty:
        return pd.DataFrame(
            {
                "distance_to_tumour_boundary_um": np.full(n, np.nan),
                "signed_distance_to_tumour_boundary_um": np.full(n, np.nan),
                "is_inside_tumour_region": np.zeros(n, dtype=bool),
            }
        )

    points = shapely.points(x, y)
    boundary = tumour_mask_polygon.boundary
    distance = shapely.distance(points, boundary)
    is_inside = shapely.contains(tumour_mask_polygon, points)
    signed_distance = np.where(is_inside, -distance, distance)

    return pd.DataFrame(
        {
            "distance_to_tumour_boundary_um": distance,
            "signed_distance_to_tumour_boundary_um": signed_distance,
            "is_inside_tumour_region": is_inside,
        }
    )


def build_tumour_boundaries(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    tumour_masks_dir = project_root / "data" / "derived" / "tumour_masks"
    output_path = project_root / "data" / "derived" / "tumour_boundaries.parquet"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )
    if not tumour_masks_dir.is_dir():
        raise PipelineError(
            f"'{tumour_masks_dir}' not found. Run `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py` first."
        )

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    obs = adata.obs[["section_id", "x_centroid", "y_centroid"]]

    section_results = []
    section_rows = []
    for section_id, group in obs.groupby("section_id", observed=True):
        mask_path = tumour_masks_dir / f"{section_id}.parquet"
        if not mask_path.is_file():
            raise PipelineError(
                f"'{mask_path}' not found. Run `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py` first."
            )
        mask_df = pd.read_parquet(mask_path)

        in_region = mask_df[mask_df["in_tumour_region"]]
        in_region_coords = obs.loc[obs.index.intersection(in_region.index)]
        polygon = build_tumour_mask_polygon(
            in_region_coords["x_centroid"].to_numpy(), in_region_coords["y_centroid"].to_numpy()
        )

        distances = compute_signed_distances(
            group["x_centroid"].to_numpy(), group["y_centroid"].to_numpy(), polygon
        )
        distances.index = group.index
        section_results.append(distances)

        section_rows.append(
            {
                "section_id": section_id,
                "n_cells": len(group),
                "n_inside_tumour_region": int(distances["is_inside_tumour_region"].sum()),
                "has_tumour_mask": polygon is not None and not polygon.is_empty,
                "median_distance_um": float(distances["distance_to_tumour_boundary_um"].median()),
            }
        )

    result = pd.concat(section_results)
    result = result.reindex(adata.obs_names)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    summary_df = pd.DataFrame(section_rows)
    return {
        "n_cells": len(result),
        "n_sections": len(summary_df),
        "n_sections_with_tumour_mask": int(summary_df["has_tumour_mask"].sum()),
        "n_cells_inside_tumour_region_total": int(summary_df["n_inside_tumour_region"].sum()),
        "fraction_cells_inside_tumour_region": round(
            float(summary_df["n_inside_tumour_region"].sum() / summary_df["n_cells"].sum()), 4
        ),
        "median_distance_by_section": summary_df.set_index("section_id")["median_distance_um"]
        .round(2)
        .to_dict(),
        "output_path": str(output_path),
    }
