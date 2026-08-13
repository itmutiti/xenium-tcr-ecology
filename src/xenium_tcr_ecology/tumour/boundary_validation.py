"""Tumour boundary validation against morphology (`07_tumour_epithelium_characterisation/06_validate_boundaries_against_morphology.py`).

Same data-availability constraint already established for
`07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py` (`morphology_concordance.py`'s module docstring) and `06_cell_type_annotation/07_blinded_annotation_review.py`'s blinded
review: this dataset holds no pathologist-drawn tumour/normal annotation,
only raw DAPI/boundary-stain morphology images, so a quantitative
"manual-review agreement" figure cannot be computed and is not
fabricated here. What is built: panels showing `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s actual
tumour-boundary polygon outline overlaid on the DAPI image, centred
on points sampled along the boundary itself (not arbitrary tissue crops --
the boundary line is where morphological agreement is actually
informative to check), plus a correctly-structured, empty review log
template for a human to fill in. That review is a human-in-the-
loop step, not performed here, for the same reason `06_cell_type_annotation/07_blinded_annotation_review.py`'s
adjudication log was not pre-filled.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd
from matplotlib.backends.backend_pdf import PdfPages

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.tumour_boundaries import build_tumour_mask_polygon

PANEL_HALF_WIDTH_UM = 75.0
N_BOUNDARY_POINTS_PER_SECTION = 4
REVIEW_LOG_COLUMNS = [
    "panel_id",
    "reviewer_id",
    "reviewer_agrees_with_boundary",
    "reviewer_confidence",
    "reviewer_notes",
    "review_date",
]


def anonymize_panel_id(
    section_id: str,
    index: int,
    salt: str = "`07_tumour_epithelium_characterisation/06_validate_boundaries_against_morphology.py`-boundary-review",
) -> str:
    return hashlib.sha256(f"{salt}:{section_id}:{index}".encode()).hexdigest()[:12]


def sample_boundary_points(
    polygon, n_points: int = N_BOUNDARY_POINTS_PER_SECTION
) -> list[tuple[float, float]]:
    """Real points sampled at even arc-length intervals ALONG the tumour
    mask's boundary line -- the boundary itself is where a morphology
    cross-check is actually informative, not an arbitrary interior or
    background crop."""
    boundary = polygon.boundary
    length = boundary.length
    if length == 0:
        return []
    points = []
    for i in range(n_points):
        frac = (i + 0.5) / n_points
        point = boundary.interpolate(frac * length)
        points.append((point.x, point.y))
    return points


def render_boundary_panel(
    sdata: sd.SpatialData,
    x: float,
    y: float,
    polygon,
    output_ax,
    half_width_um: float = PANEL_HALF_WIDTH_UM,
) -> None:
    min_coord = [x - half_width_um, y - half_width_um]
    max_coord = [x + half_width_um, y + half_width_um]

    cropped = sd.bounding_box_query(
        sdata,
        axes=("x", "y"),
        min_coordinate=min_coord,
        max_coordinate=max_coord,
        target_coordinate_system="global",
    )
    if cropped is None:
        raise PipelineError(f"Bounding box query at ({x}, {y}) returned no data.")

    image_key = next(iter(cropped.images), None)
    if image_key is not None:
        img = cropped.images[image_key]
        img_data = np.asarray(img.isel(c=0) if "c" in img.dims else img)
        # Same fix as `06_cell_type_annotation/07_blinded_annotation_review.py`, `07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`: the cropped image resets to a local pixel frame,
        # while the polygon remains in physical/global coordinates.
        output_ax.imshow(
            img_data,
            cmap="gray",
            origin="lower",
            extent=[min_coord[0], max_coord[0], min_coord[1], max_coord[1]],
        )

    from shapely.geometry import box

    crop_box = box(min_coord[0], min_coord[1], max_coord[0], max_coord[1])
    visible_boundary = polygon.boundary.intersection(crop_box)
    if not visible_boundary.is_empty:
        geoms = (
            list(visible_boundary.geoms)
            if hasattr(visible_boundary, "geoms")
            else [visible_boundary]
        )
        for geom in geoms:
            if geom.geom_type == "LineString":
                xs, ys = geom.xy
                output_ax.plot(xs, ys, color="cyan", linewidth=1.5)

    output_ax.set_xlim(min_coord[0], max_coord[0])
    output_ax.set_ylim(min_coord[1], max_coord[1])
    output_ax.set_xticks([])
    output_ax.set_yticks([])


def build_boundary_validation_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    tumour_masks_dir = project_root / "data" / "derived" / "tumour_masks"
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_path = project_root / "reports" / "tumour" / "boundary_validation.pdf"
    review_log_path = project_root / "reports" / "tumour" / "boundary_review_log.tsv"

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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_rendered = 0
    n_failed = 0
    n_sections_no_boundary = 0
    log_rows = []

    with PdfPages(output_path) as pdf:
        for section_id, group in obs.groupby("section_id", observed=True):
            mask_path = tumour_masks_dir / f"{section_id}.parquet"
            if not mask_path.is_file():
                raise PipelineError(
                    f"'{mask_path}' not found. Run `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py` first."
                )
            mask_df = pd.read_parquet(mask_path)
            in_region = mask_df[mask_df["in_tumour_region"]]
            in_region_coords = group.loc[group.index.intersection(in_region.index)]
            polygon = build_tumour_mask_polygon(
                in_region_coords["x_centroid"].to_numpy(), in_region_coords["y_centroid"].to_numpy()
            )
            if polygon is None or polygon.is_empty:
                n_sections_no_boundary += 1
                continue

            boundary_points = sample_boundary_points(polygon)
            zarr_path = spatialdata_root / f"{section_id}.zarr"
            if not zarr_path.exists():
                n_failed += len(boundary_points)
                continue
            sdata = sd.read_zarr(zarr_path)

            for i, (x, y) in enumerate(boundary_points):
                panel_id = anonymize_panel_id(section_id, i)
                fig, ax = plt.subplots(figsize=(5, 5))
                try:
                    render_boundary_panel(sdata, x, y, polygon, ax)
                    ax.set_title(panel_id, fontsize=13, fontfamily="Liberation Sans")
                    n_rendered += 1
                except PipelineError:
                    n_failed += 1
                    plt.close(fig)
                    continue
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                log_rows.append({"panel_id": panel_id, "section_id": section_id, "x": x, "y": y})

    log_df = pd.DataFrame(log_rows)
    for col in REVIEW_LOG_COLUMNS:
        if col not in log_df.columns and col != "panel_id":
            log_df[col] = pd.NA
    review_log_path.parent.mkdir(parents=True, exist_ok=True)
    log_df[["panel_id"] + [c for c in REVIEW_LOG_COLUMNS if c != "panel_id"]].to_csv(
        review_log_path, sep="\t", index=False
    )

    return {
        "n_panels_rendered": n_rendered,
        "n_panels_failed": n_failed,
        "n_sections_no_boundary": n_sections_no_boundary,
        "output_path": str(output_path),
        "review_log_path": str(review_log_path),
        "review_status": f"PENDING_HUMAN_REVIEW -- 0 of {n_rendered} panel(s) reviewed",
    }
