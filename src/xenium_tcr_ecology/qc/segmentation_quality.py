"""Segmentation quality assessment (`04_quality_control/03_assess_segmentation_quality.py`).

Real geometric checks against the actual cell/nucleus boundary polygons
(not just the pre-computed area columns from cells.parquet, which cannot
detect e.g. a nucleus polygon that extends outside its cell polygon):

  - polygon validity (shapely .is_valid -- self-intersecting or otherwise
    malformed polygons indicate a segmentation/export error)
  - nucleus-within-cell overlap (biologically, a nucleus should sit inside
    its own cell). Measured as a *fractional* area overlap
    (intersection(cell, nucleus).area / nucleus.area), not strict polygon
    containment: Xenium's cell boundaries are constructed by expanding
    outward from the nucleus boundary, so the two polygons legitimately
    share edges/vertices. Empirically (checked against the data, n=1000
    sampled cells), shapely's strict `.contains()` flags ~43% of cells as
    "not contained" purely from sub-pixel floating-point slivers poking past
    a shared boundary, even though the true overlap fraction for those same
    cells is >=0.974. A fractional-overlap threshold (see
    NUCLEUS_OVERLAP_FLAG_THRESHOLD) is therefore the correct check; the
    threshold sits comfortably below the observed floating-point-artifact
    floor (0.974) so it does not re-introduce false flags, while still
    catching genuine segmentation errors (nucleus meaningfully outside its
    cell).
  - nucleus:cell area ratio > 1 (geometrically impossible if segmentation
    is correct, since a nucleus cannot be larger than the cell containing
    it -- a clear, unambiguous error indicator)

Runs on a sample of cells per section (shapely polygon operations on all
~1.19M cells would be slow and are not necessary to estimate a section's
segmentation quality) -- matches the blueprint's own "samples cells for
visual review" framing for this phase.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed

DEFAULT_SAMPLE_SIZE = 1000

# Below the empirically observed floating-point-artifact floor of 0.974 (see
# module docstring) -- flags genuine containment failures without
# re-triggering on shared-boundary sub-pixel slivers.
NUCLEUS_OVERLAP_FLAG_THRESHOLD = 0.95


def check_cell_nucleus_pair(cell_poly, nucleus_poly) -> dict:
    """Pure geometric check for one cell/nucleus polygon pair -- factored out
    so the core logic is testable with plain shapely polygons, not only via
    a full SpatialData store."""
    cell_valid = cell_poly.is_valid
    nucleus_valid = nucleus_poly.is_valid
    overlap_fraction = None
    contained = False
    area_exceeds = False
    if cell_valid and nucleus_valid and nucleus_poly.area > 0:
        overlap_fraction = cell_poly.intersection(nucleus_poly).area / nucleus_poly.area
        contained = overlap_fraction >= NUCLEUS_OVERLAP_FLAG_THRESHOLD
        area_exceeds = nucleus_poly.area > cell_poly.area
    return {
        "cell_valid": cell_valid,
        "nucleus_valid": nucleus_valid,
        "nucleus_overlap_fraction": overlap_fraction,
        "nucleus_contained": contained,
        "nucleus_area_exceeds_cell": area_exceeds,
    }


def assess_section_segmentation(
    zarr_path: Path, sample_size: int = DEFAULT_SAMPLE_SIZE, rng_seed: int = get_default_seed()
) -> dict:
    if not zarr_path.exists():
        raise PipelineError(
            f"'{zarr_path}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    sdata = sd.read_zarr(zarr_path)
    cell_shapes = sdata["cell_boundaries"]
    nucleus_shapes = sdata["nucleus_boundaries"]

    common_ids = cell_shapes.index.intersection(nucleus_shapes.index)
    if len(common_ids) == 0:
        raise PipelineError(f"'{zarr_path}': no cells with both cell and nucleus boundaries.")

    rng = np.random.default_rng(rng_seed)
    sample_ids = rng.choice(common_ids, size=min(sample_size, len(common_ids)), replace=False)

    n_invalid_cell_polygon = 0
    n_invalid_nucleus_polygon = 0
    n_nucleus_not_contained = 0
    n_nucleus_area_exceeds_cell = 0
    overlap_fractions: list[float] = []
    example_flagged: list[str] = []

    for cid in sample_ids:
        cell_poly = cell_shapes.loc[cid, "geometry"]
        nucleus_poly = nucleus_shapes.loc[cid, "geometry"]
        result = check_cell_nucleus_pair(cell_poly, nucleus_poly)

        if not result["cell_valid"]:
            n_invalid_cell_polygon += 1
        if not result["nucleus_valid"]:
            n_invalid_nucleus_polygon += 1

        flagged = False
        if result["cell_valid"] and result["nucleus_valid"]:
            overlap_fractions.append(result["nucleus_overlap_fraction"])
            if not result["nucleus_contained"]:
                n_nucleus_not_contained += 1
                flagged = True
            if result["nucleus_area_exceeds_cell"]:
                n_nucleus_area_exceeds_cell += 1
                flagged = True
        if flagged and len(example_flagged) < 10:
            example_flagged.append(str(cid))

    n = len(sample_ids)
    return {
        "n_cells_sampled": n,
        "fraction_invalid_cell_polygon": round(n_invalid_cell_polygon / n, 4),
        "fraction_invalid_nucleus_polygon": round(n_invalid_nucleus_polygon / n, 4),
        "median_nucleus_overlap_fraction": (
            round(float(np.median(overlap_fractions)), 4) if overlap_fractions else None
        ),
        "fraction_nucleus_not_contained": round(n_nucleus_not_contained / n, 4),
        "fraction_nucleus_area_exceeds_cell": round(n_nucleus_area_exceeds_cell / n, 4),
        "example_flagged_cell_ids": ";".join(example_flagged),
    }


def build_segmentation_quality_report(spatialdata_root: Path, output_path: Path) -> dict:
    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        raise PipelineError(
            f"No .zarr stores found under '{spatialdata_root}'. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    rows = []
    for zarr_path in zarr_paths:
        section_id = zarr_path.stem
        metrics = assess_section_segmentation(zarr_path)
        rows.append({"section_id": section_id, **metrics})

    df = pd.DataFrame(rows)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Machine-readable per-section summary alongside the HTML review, so
    # downstream consumers (`04_quality_control/09_generate_qc_release_report.py`'s QC release report) do not need to
    # re-run this section's ~1000-cells-per-section shapely computation or
    # scrape the HTML table just to read back numbers already computed here.
    df.drop(columns=["example_flagged_cell_ids"]).to_parquet(output_path.with_suffix(".parquet"))
    rows_html = "".join(
        f"<tr><td>{r['section_id']}</td><td>{r['n_cells_sampled']}</td>"
        f"<td>{r['fraction_invalid_cell_polygon']:.4f}</td>"
        f"<td>{r['fraction_invalid_nucleus_polygon']:.4f}</td>"
        f"<td>{r['median_nucleus_overlap_fraction']:.4f}</td>"
        f"<td>{r['fraction_nucleus_not_contained']:.4f}</td>"
        f"<td>{r['fraction_nucleus_area_exceeds_cell']:.4f}</td>"
        f"<td>{r['example_flagged_cell_ids']}</td></tr>"
        for r in rows
    )
    html = (
        "<html><head><title>Segmentation quality review</title></head><body>"
        f"<p>Sampled-cell geometric QC: polygon validity, nucleus-within-cell area overlap "
        f"(fractional, threshold {NUCLEUS_OVERLAP_FLAG_THRESHOLD} -- not strict polygon "
        f"containment, since Xenium cell/nucleus boundaries legitimately share edges "
        f"from the nucleus-expansion construction), nucleus:cell area ratio. Flags "
        f"candidates for review; does not exclude anything (`04_quality_control/06_define_qc_thresholds_hierarchically.R`, `04_quality_control/07_apply_qc_filters_with_audit_trail.py`'s job).</p>"
        "<table border='1' cellpadding='6' cellspacing='0'>"
        "<tr><th>Section</th><th>N sampled</th><th>Frac. invalid cell polygon</th>"
        "<th>Frac. invalid nucleus polygon</th><th>Median nucleus overlap fraction</th>"
        "<th>Frac. nucleus not contained</th><th>Frac. nucleus area &gt; cell area</th>"
        "<th>Example flagged cell IDs</th></tr>"
        f"{rows_html}</table></body></html>"
    )
    output_path.write_text(html)

    has_invalid_polygons = (df["fraction_invalid_cell_polygon"] > 0) | (
        df["fraction_invalid_nucleus_polygon"] > 0
    )
    return {
        "sections_processed": len(zarr_paths),
        "median_fraction_nucleus_not_contained": float(
            df["fraction_nucleus_not_contained"].median()
        ),
        "max_fraction_nucleus_not_contained": float(df["fraction_nucleus_not_contained"].max()),
        "sections_with_invalid_polygons": int(has_invalid_polygons.sum()),
    }
