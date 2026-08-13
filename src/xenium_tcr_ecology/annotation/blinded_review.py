"""Blinded spatial review panel generation (`06_cell_type_annotation/07_blinded_annotation_review.py`).

Generates blinded (predicted-label-hidden) spatial context panels for
expert review, stratified across sections, major lineages, and confidence
tiers, deliberately oversampling low-confidence/ambiguous cells (Phase
6.06) since those most need expert scrutiny.

This module builds everything that can be built without a human domain
expert: the panels themselves (a spatial crop of the DAPI morphology
image with cell/nucleus boundary outlines, no colour-coding by predicted
type) and a correctly-structured, empty adjudication log template. It does
not fabricate adjudication decisions -- judging whether a specific cell's
tissue morphology matches its algorithmic prediction requires a
pathologist/immunologist actually looking at each panel, which this
pipeline cannot supply. Completing `06_cell_type_annotation/07_blinded_annotation_review.py` (i.e. filling in
`adjudication_log.tsv`'s adjudication columns) is therefore a genuine,
explicit human-in-the-loop step, not a further automation gap.

Panel IDs are an anonymised hash of the cell ID, not the cell ID itself
(blinding): the mapping back to the cell ID and its algorithmic
prediction is kept in a separate key file, not shown to the reviewer.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError

PANEL_HALF_WIDTH_UM = 50.0
N_REVIEW_PANELS = 150
# A documented judgment call: ambiguous cells are the ones most in need
# of expert scrutiny, so they are deliberately over-represented in the
# review sample relative to their 25.3% population share (`06_cell_type_annotation/06_integrate_annotation_evidence.py`).
AMBIGUOUS_OVERSAMPLE_TARGET_FRACTION = 0.5
ADJUDICATION_LOG_COLUMNS = [
    "panel_id",
    "reviewer_id",
    "reviewer_adjudicated_lineage",
    "reviewer_confidence_in_own_call",
    "reviewer_notes",
    "adjudication_date",
]


def anonymize_panel_id(
    cell_id: str,
    salt: str = "`06_cell_type_annotation/07_blinded_annotation_review.py`-blinded-review",
) -> str:
    return hashlib.sha256(f"{salt}:{cell_id}".encode()).hexdigest()[:12]


def select_review_sample(
    final_annotations: pd.DataFrame,
    n_total: int = N_REVIEW_PANELS,
    rng_seed: int = get_annotation_seed(),
) -> pd.DataFrame:
    if len(final_annotations) == 0:
        raise PipelineError("final_annotations is empty -- nothing to sample.")

    rng = np.random.default_rng(rng_seed)
    n_ambiguous_target = int(round(n_total * AMBIGUOUS_OVERSAMPLE_TARGET_FRACTION))
    n_confident_target = n_total - n_ambiguous_target

    def _stratified_sample(pool: pd.DataFrame, n_target: int) -> pd.DataFrame:
        if len(pool) == 0:
            return pool.iloc[0:0]
        lineages = pool["final_lineage"].unique()
        per_lineage = max(1, n_target // len(lineages))
        parts = []
        for lineage in lineages:
            sub = pool[pool["final_lineage"] == lineage]
            n = min(per_lineage, len(sub))
            idx = rng.choice(sub.index, size=n, replace=False)
            parts.append(sub.loc[idx])
        return pd.concat(parts)

    ambiguous_sample = _stratified_sample(
        final_annotations[final_annotations["is_ambiguous"]], n_ambiguous_target
    )
    confident_sample = _stratified_sample(
        final_annotations[~final_annotations["is_ambiguous"]], n_confident_target
    )
    return pd.concat([ambiguous_sample, confident_sample])


def render_panel(
    sdata: sd.SpatialData,
    x: float,
    y: float,
    output_path: Path,
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

    fig, ax = plt.subplots(figsize=(5, 5))
    image_key = next(iter(cropped.images), None)
    if image_key is not None:
        img = cropped.images[image_key]
        img_data = np.asarray(img.isel(c=0) if "c" in img.dims else img)
        # sd.bounding_box_query() resets a cropped image's own x/y
        # coordinates to a local pixel-index frame (confirmed:
        # 0.5-470.5 for a 100um crop), while shapes remain in the original
        # physical/global coordinate system -- plotting both without an
        # explicit, physically-correct `extent` produced a bug on the
        # first render (the image collapsed into a tiny corner, disconnected
        # from the cell boundary outlines). The requested bounding box
        # itself (min_coord/max_coord, already in the same physical/global
        # system as the shapes) is the correct extent, not the image's own
        # post-crop coordinate labels.
        ax.imshow(
            img_data,
            cmap="gray",
            origin="lower",
            extent=[min_coord[0], max_coord[0], min_coord[1], max_coord[1]],
        )

    # Cell/nucleus boundaries drawn as neutral outlines only -- no
    # colour-coding by predicted lineage, which would break blinding.
    for shapes_key, color in [("cell_boundaries", "cyan"), ("nucleus_boundaries", "yellow")]:
        if shapes_key in cropped.shapes:
            cropped.shapes[shapes_key].boundary.plot(ax=ax, color=color, linewidth=0.5)

    ax.set_xlim(min_coord[0], max_coord[0])
    ax.set_ylim(min_coord[1], max_coord[1])
    ax.set_title("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def build_blinded_review_panels(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    panels_dir = project_root / "reports" / "annotation" / "blinded_panels"
    key_path = project_root / "reports" / "annotation" / "blinded_panel_key.tsv"
    adjudication_log_path = project_root / "reports" / "annotation" / "adjudication_log.tsv"

    if not matrix_path.is_file():
        raise PipelineError(
            f"'{matrix_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )
    if not final_annotations_path.is_file():
        raise PipelineError(
            f"'{final_annotations_path}' not found. Run `06_cell_type_annotation/06_integrate_annotation_evidence.py` first."
        )

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    final_annotations = pd.read_parquet(final_annotations_path)
    final_annotations = final_annotations.join(
        adata.obs[["section_id", "x_centroid", "y_centroid"]]
    )

    sample = select_review_sample(final_annotations)

    panels_dir.mkdir(parents=True, exist_ok=True)
    key_rows = []
    n_rendered = 0
    n_failed = 0
    for section_id, group in sample.groupby("section_id", observed=True):
        zarr_path = spatialdata_root / f"{section_id}.zarr"
        if not zarr_path.exists():
            n_failed += len(group)
            continue
        sdata = sd.read_zarr(zarr_path)
        for cell_id, row in group.iterrows():
            panel_id = anonymize_panel_id(cell_id)
            output_path = panels_dir / f"{panel_id}.png"
            try:
                render_panel(sdata, row["x_centroid"], row["y_centroid"], output_path)
                n_rendered += 1
            except PipelineError:
                n_failed += 1
                continue
            key_rows.append(
                {
                    "panel_id": panel_id,
                    "cell_id": cell_id,
                    "section_id": section_id,
                    "algorithmic_final_lineage": row["final_lineage"],
                    "algorithmic_final_substate": row["final_substate"],
                    "algorithmic_confidence": row["confidence"],
                    "algorithmic_is_ambiguous": row["is_ambiguous"],
                }
            )

    key_df = pd.DataFrame(key_rows)
    key_df.to_csv(key_path, sep="\t", index=False)

    # Empty template, correct schema, ready for a human reviewer -- no
    # adjudication rows are fabricated here.
    adjudication_log = pd.DataFrame(columns=ADJUDICATION_LOG_COLUMNS)
    adjudication_log["panel_id"] = key_df["panel_id"]
    adjudication_log.to_csv(adjudication_log_path, sep="\t", index=False)

    return {
        "n_panels_requested": N_REVIEW_PANELS,
        "n_panels_rendered": n_rendered,
        "n_panels_failed": n_failed,
        "n_ambiguous_in_sample": int(sample["is_ambiguous"].sum()),
        "panels_dir": str(panels_dir),
        "key_path": str(key_path),
        "adjudication_log_path": str(adjudication_log_path),
        "adjudication_status": "PENDING_HUMAN_REVIEW -- 0 of {} panels adjudicated".format(
            n_rendered
        ),
    }
