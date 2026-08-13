"""Spatial compartment labelling: tumour core, inner margin, outer margin,
distal stroma (`07_tumour_epithelium_characterisation/07_define_invasive_front_and_compartments.py`).

Labels every cell (all lineages, matching `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s "for all
cells" scope) using `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s `signed_distance_to_tumour_boundary_um`
and a predeclared distance band, per the specification's "predeclared
distance bands and sensitivity analyses" requirement.

**Band width, a judgment call grounded directly in the achieved data,
not imported unchanged from the literature:** checked the distribution
before choosing anything. Common spatial-oncology invasive-margin
definitions (e.g. the barrier-topology literature already cited in this
project's registry, governance/analysis_registry.tsv's
`q3_literature_benchmark`) use 25-100um bands -- but at every one of
those widths, this dataset's `tumour_core` fraction is exactly 0%
(confirmed: the maximum achieved interior depth,
`signed_distance_to_tumour_boundary_um.min()`, is only -14.55um across
all 1,094,816 non-null cells). This is a structural consequence of
`07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py`, `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s boundary construction: the tumour mask is a buffered
union of individual cell points (`BOUNDARY_BUFFER_UM = 5.0`), so its
boundary hugs individual cell clusters rather than approximating one
smoothed macro-scale tumour margin -- there is very little "deep interior"
space for a classical 25-100um core zone to exist in, by construction,
not because tumour cells are absent. `PRIMARY_BAND_WIDTH_UM = 3.0` was
chosen instead because it is the smallest, still-round band width that
resolves a non-degenerate `tumour_core` population under the
achievable depth range (checked: at band=3um, 7.27% of cells
are `tumour_core`; at band=5um, only 0.02% are, because the bulk of
interior points sit within the ~5um buffer radius itself, not beyond
it).

**Resolution caveat:** `tumour_core`/`inner_margin` at
this resolution reflect fine local geometry (how close a cell is to its
own small cluster's edge) rather than a classical "deep tumour mass core",
an inherited limitation from the per-cell (not macro-region) boundary
construction, not a new one introduced here
(see e.g. `06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s squamous-marker-coverage limitation).

**Sensitivity analysis, per the specification's explicit requirement:**
compartments are also computed at `SENSITIVITY_BAND_WIDTHS_UM = [1.5,
6.0]` (half and double the primary), and all three assignments are
retained in the output so downstream consumers can check whether a
finding is band-width-sensitive, not only reported once.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

PRIMARY_BAND_WIDTH_UM = 3.0
SENSITIVITY_BAND_WIDTHS_UM = [1.5, 6.0]

COMPARTMENT_ORDER = ["tumour_core", "inner_margin", "outer_margin", "distal_stroma"]


def assign_compartment(signed_distance: pd.Series, band_width_um: float) -> pd.Series:
    """Pure, testable compartment assignment for one band width. NaN
    signed_distance (no tumour mask in that section, `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`) maps to
    NaN compartment -- a "no tumour boundary in this section" result,
    not imputed or dropped."""
    conditions = [
        signed_distance <= -band_width_um,
        (signed_distance > -band_width_um) & (signed_distance <= 0),
        (signed_distance > 0) & (signed_distance <= band_width_um),
        signed_distance > band_width_um,
    ]
    result = pd.Series(
        np.select(conditions, COMPARTMENT_ORDER, default=None), index=signed_distance.index
    )
    result[signed_distance.isna()] = np.nan
    return result


def build_spatial_compartments(project_root: Path) -> dict:
    tumour_boundaries_path = project_root / "data" / "derived" / "tumour_boundaries.parquet"
    output_path = project_root / "data" / "derived" / "spatial_compartments.parquet"

    if not tumour_boundaries_path.is_file():
        raise PipelineError(
            f"'{tumour_boundaries_path}' not found. Run `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py` first."
        )

    boundaries = pd.read_parquet(tumour_boundaries_path)
    signed_distance = boundaries["signed_distance_to_tumour_boundary_um"]

    result = pd.DataFrame(index=boundaries.index)
    result["compartment"] = assign_compartment(signed_distance, PRIMARY_BAND_WIDTH_UM)
    for band in SENSITIVITY_BAND_WIDTHS_UM:
        result[f"compartment_band{band}um"] = assign_compartment(signed_distance, band)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    primary_counts = result["compartment"].value_counts()
    n_valid = result["compartment"].notna().sum()
    sensitivity_fractions = {}
    for band in [PRIMARY_BAND_WIDTH_UM] + SENSITIVITY_BAND_WIDTHS_UM:
        col = "compartment" if band == PRIMARY_BAND_WIDTH_UM else f"compartment_band{band}um"
        counts = result[col].value_counts()
        n = result[col].notna().sum()
        sensitivity_fractions[f"band_{band}um"] = {
            c: round(float(counts.get(c, 0) / n), 4) for c in COMPARTMENT_ORDER
        }

    return {
        "n_cells": len(result),
        "n_cells_with_compartment": int(n_valid),
        "n_cells_no_tumour_mask": int(result["compartment"].isna().sum()),
        "primary_band_width_um": PRIMARY_BAND_WIDTH_UM,
        "sensitivity_band_widths_um": SENSITIVITY_BAND_WIDTHS_UM,
        "primary_compartment_counts": {k: int(v) for k, v in primary_counts.items()},
        "compartment_fractions_by_band_width": sensitivity_fractions,
        "output_path": str(output_path),
    }
