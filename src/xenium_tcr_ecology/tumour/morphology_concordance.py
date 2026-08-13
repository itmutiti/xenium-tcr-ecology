"""Malignancy-call cross-validation against morphology (`07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py`).

The specification calls for cross-checking the transcriptional
malignancy call against "any co-registered morphology/pathology
annotation." Checked before implementing: this
project's data holdings (`data/standardised/*/morphology.ome.tif.gz`,
`data/staged/*/*_morphology.ome.tif.gz`) contain only raw DAPI/boundary-
stain morphology images -- there is no pathologist-drawn tumour/normal
region annotation anywhere in this dataset (confirmed by a search
of `data/` and `config/metadata/` before writing this module). A
quantitative sensitivity/specificity concordance against a ground truth
therefore cannot be computed here, and is not fabricated.

This is the same human-in-the-loop data-availability gap already
encountered once in this project (`06_cell_type_annotation/07_blinded_annotation_review.py`'s blinded review):
what can be built without a
domain expert is built; what requires one is not simulated. Concretely,
this module produces two non-fabricated things instead:

1. **A quantitative spatial-coherence sanity check** (`compute_malignancy_
   spatial_autocorrelation`): per-section Moran's I of `malignancy_
   probability` over a spatial k-NN graph (squidpy). This is a genuine,
   partial circularity-risk bound, not a full substitute for
   pathologist validation: a malignancy score dominated by per-cell
   noise/technical artefacts (the kind of thing an uninformative or
   circular scoring scheme would produce) would show weak spatial
   autocorrelation, since noise has no reason to be spatially clustered;
   tumour tissue is spatially contiguous, so a genuine malignancy signal
   should show strong positive autocorrelation. A high Moran's I is
   evidence the score reflects tissue-level structure; it is evidence
   against pure noise/circularity, not proof the structure it reflects
   is correctly labelled "malignant" specifically.
2. **Visual overlay panels** (`render_malignancy_overlay_panel`),
   reusing `06_cell_type_annotation/07_blinded_annotation_review.py`'s established DAPI-crop rendering fix (the
   `sd.bounding_box_query()` local-pixel-vs-physical-coordinate mismatch): cells coloured
   by `malignancy_probability`, overlaid on the DAPI image, for a sample
   of high- and low-malignancy-probability regions per patient -- for a
   human to actually look at and judge whether high-probability regions
   show morphology consistent with malignancy (nuclear atypia,
   disorganised architecture). That judgement itself is not performed
   here, for the same reason `06_cell_type_annotation/07_blinded_annotation_review.py`'s adjudication log is not
   pre-filled: it requires a human domain expert.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd
import squidpy as sq
from matplotlib.backends.backend_pdf import PdfPages

from xenium_tcr_ecology.infra.exceptions import PipelineError

SPATIAL_KNN_NEIGHBORS = 6
PANEL_HALF_WIDTH_UM = 75.0
N_EXTREME_REGIONS_PER_PATIENT = 2


def compute_malignancy_spatial_autocorrelation(
    x: np.ndarray,
    y: np.ndarray,
    malignancy_probability: np.ndarray,
    n_neighs: int = SPATIAL_KNN_NEIGHBORS,
) -> dict:
    """Pure, testable Moran's I computation for one section's worth of
    cells -- factored out so it is testable with plain coordinate arrays."""
    n = len(x)
    if n < n_neighs + 1:
        return {"n_cells": n, "morans_i": None, "pval_norm": None}

    sub = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32))
    sub.obsm["spatial"] = np.column_stack([x, y])
    sub.obs["malignancy_probability"] = malignancy_probability

    sq.gr.spatial_neighbors_knn(sub, n_neighs=n_neighs)
    result = sq.gr.spatial_autocorr(
        sub,
        mode="moran",
        attr="obs",
        genes=["malignancy_probability"],
        copy=True,
        show_progress_bar=False,
    )
    row = result.loc["malignancy_probability"]
    return {"n_cells": n, "morans_i": float(row["I"]), "pval_norm": float(row["pval_norm"])}


def render_malignancy_overlay_panel(
    sdata: sd.SpatialData,
    x: float,
    y: float,
    cell_coords: np.ndarray,
    malignancy_values: np.ndarray,
    ax,
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
        # Same fix as `06_cell_type_annotation/07_blinded_annotation_review.py`: the cropped image's own coordinates reset to a local
        # pixel frame, while cell coordinates remain physical/global -- the
        # originally-requested bounding box is the correct plotting extent.
        ax.imshow(
            img_data,
            cmap="gray",
            origin="lower",
            extent=[min_coord[0], max_coord[0], min_coord[1], max_coord[1]],
        )

    in_view = (
        (cell_coords[:, 0] >= min_coord[0])
        & (cell_coords[:, 0] <= max_coord[0])
        & (cell_coords[:, 1] >= min_coord[1])
        & (cell_coords[:, 1] <= max_coord[1])
    )
    sc = ax.scatter(
        cell_coords[in_view, 0],
        cell_coords[in_view, 1],
        c=malignancy_values[in_view],
        cmap="viridis",
        vmin=0,
        vmax=1,
        s=4,
        alpha=0.8,
    )
    ax.set_xlim(min_coord[0], max_coord[0])
    ax.set_ylim(min_coord[1], max_coord[1])
    ax.set_xticks([])
    ax.set_yticks([])
    return sc


def select_extreme_regions(
    section_scores: pd.DataFrame, n_per_side: int = N_EXTREME_REGIONS_PER_PATIENT
) -> pd.DataFrame:
    """One highest- and one lowest-malignancy_probability cell per patient
    section (used as panel centres) -- a deterministic, non-random
    selection of the most informative regions to visually inspect."""
    parts = []
    for section_id, group in section_scores.groupby("section_id", observed=True):
        top = group.nlargest(n_per_side, "malignancy_probability").assign(region_type="high")
        bottom = group.nsmallest(n_per_side, "malignancy_probability").assign(region_type="low")
        parts.append(pd.concat([top, bottom]))
    return pd.concat(parts)


def build_morphology_concordance_report(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    malignancy_scores_path = project_root / "data" / "derived" / "malignancy_scores.parquet"
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_path = project_root / "reports" / "tumour" / "morphology_concordance.pdf"
    autocorr_path = project_root / "reports" / "tumour" / "malignancy_spatial_autocorrelation.tsv"

    for p in (matrix_path, malignancy_scores_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    adata = ad.read_h5ad(matrix_path)
    scores = pd.read_parquet(malignancy_scores_path)
    scores = scores.join(adata.obs[["x_centroid", "y_centroid"]])

    autocorr_rows = []
    for section_id, group in scores.groupby("section_id", observed=True):
        result = compute_malignancy_spatial_autocorrelation(
            group["x_centroid"].to_numpy(),
            group["y_centroid"].to_numpy(),
            group["malignancy_probability"].to_numpy(),
        )
        autocorr_rows.append({"section_id": section_id, **result})
    autocorr_df = pd.DataFrame(autocorr_rows)
    autocorr_path.parent.mkdir(parents=True, exist_ok=True)
    autocorr_df.to_csv(autocorr_path, sep="\t", index=False)

    regions = select_extreme_regions(scores)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_rendered = 0
    n_failed = 0
    with PdfPages(output_path) as pdf:
        for section_id, section_regions in regions.groupby("section_id", observed=True):
            zarr_path = spatialdata_root / f"{section_id}.zarr"
            if not zarr_path.exists():
                n_failed += len(section_regions)
                continue
            sdata = sd.read_zarr(zarr_path)
            section_cells = scores[scores["section_id"] == section_id]
            cell_coords = section_cells[["x_centroid", "y_centroid"]].to_numpy()
            malignancy_values = section_cells["malignancy_probability"].to_numpy()

            for _, row in section_regions.iterrows():
                fig, ax = plt.subplots(figsize=(5, 5))
                try:
                    scatter = render_malignancy_overlay_panel(
                        sdata,
                        row["x_centroid"],
                        row["y_centroid"],
                        cell_coords,
                        malignancy_values,
                        ax,
                    )
                    ax.set_title(
                        f"{section_id} -- {row['region_type']} malignancy region",
                        fontsize=13,
                        fontfamily="Liberation Sans",
                    )
                    if scatter is not None:
                        fig.colorbar(scatter, ax=ax, label="Malignancy probability", fraction=0.046)
                    n_rendered += 1
                except PipelineError:
                    n_failed += 1
                    plt.close(fig)
                    continue
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

    return {
        "n_sections_with_autocorrelation": int(autocorr_df["morans_i"].notna().sum()),
        "median_morans_i": round(float(autocorr_df["morans_i"].median()), 4),
        "min_morans_i": round(float(autocorr_df["morans_i"].min()), 4),
        "max_morans_i": round(float(autocorr_df["morans_i"].max()), 4),
        "n_panels_rendered": n_rendered,
        "n_panels_failed": n_failed,
        "autocorr_path": str(autocorr_path),
        "output_path": str(output_path),
        "concordance_status": (
            "Partial -- no pathologist tumour/normal annotation exists in this dataset; "
            "quantitative sensitivity/specificity cannot be computed. Spatial autocorrelation and "
            "visual overlay panels are non-fabricated partial checks; qualitative morphology judgement "
            "of the rendered panels is a human-in-the-loop step, not performed here."
        ),
    }
