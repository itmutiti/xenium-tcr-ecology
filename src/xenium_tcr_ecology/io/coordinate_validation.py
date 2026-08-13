"""Cross-section coordinate system validation (`03_spatialdata_import/02_validate_coordinate_systems.py`).

Extends the single-section physical-alignment check performed manually
during `03_spatialdata_import/01_import_each_section_to_spatialdata.py`'s development into a systematic, per-section, automated
check across every imported section.

Methodology note, established empirically while building this: comparing
cell centroids against globally-random background pixels (the first
version of this check) gave false failures on the 4 densest sections
(4800-5700 cells/mm^2) -- in dense tissue, a large fraction of "random
background" pixels land on cells too, inflating the background baseline
and compressing the apparent signal ratio even though alignment is
correct. Confirmed by checking: the 4 "failing" sections under that
method were exactly the 4 highest-density sections, with a clean monotonic
relationship between density and apparent failure -- not the pattern a
coordinate bug would produce (which has no reason to correlate with
tissue density).

Second fix attempt (single fixed shift, ~1 nuclear diameter) still gave
false failures on 3 patients (9, 12, 13, 6 sections) with a physical
explanation ruled in by direct investigation, not assumed: profiling mean
DAPI intensity across a range of shift distances (0, 2, 4, 8, 16, 32, 64 um)
for a "failing" section (P09_run1) showed a smooth, monotonically
*declining* profile starting from shift=0 as the maximum -- exactly what
correct alignment produces, just with a broader, more gradual peak than
sparser sections (consistent with these patients' tissue being unusually
densely and closely packed, not with a coordinate bug, which would show the
peak away from shift=0, or no distance-intensity relationship at all).

Final methodology: instead of one arbitrary shift distance, confirm shift=0
is the local maximum across a small profile of increasing shift distances.
This is the density-independent, physically-motivated criterion -- a
DAPI-stained nucleus peaks at its own centre regardless of how close its
neighbours are; only a genuine coordinate error would displace that peak
away from the recorded (x, y) position.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter
from xenium_tcr_ecology.infra.seeding import get_default_seed

VALIDATION_FIELDS = [
    "section_id",
    "n_cells",
    "n_cell_boundaries",
    "n_nucleus_boundaries",
    "coordinate_systems",
    "intensity_profile_by_shift_um",
    "peak_at_zero_shift",
    "status",
]

PROFILE_SHIFTS_UM = [0.0, 2.0, 4.0, 8.0, 16.0, 32.0]


def _mean_intensity_at_shift(img, sample, pixel_size_um, shift_um, rng_seed):
    rng = np.random.default_rng(rng_seed)
    shift_px = shift_um / pixel_size_um
    angles = rng.uniform(0, 2 * np.pi, size=len(sample))
    vals = []
    for (_, row), ang in zip(sample.iterrows(), angles):
        x = row["x_centroid"] / pixel_size_um + shift_px * np.cos(ang)
        y = row["y_centroid"] / pixel_size_um + shift_px * np.sin(ang)
        xi, yi = int(x), int(y)
        if 0 <= yi < img.shape[0] and 0 <= xi < img.shape[1]:
            vals.append(img[yi, xi])
    return float(np.mean(vals)) if vals else 0.0


def validate_section(
    zarr_path: Path, pixel_size_um: float, rng_seed: int = get_default_seed()
) -> dict:
    sdata = sd.read_zarr(zarr_path)

    n_cells = sdata["table"].n_obs
    n_cell_boundaries = sdata["cell_boundaries"].shape[0]
    n_nucleus_boundaries = sdata["nucleus_boundaries"].shape[0]

    img = sdata["morphology_mip"].compute().to_numpy()[0]
    sample = sdata["table"].obs.sample(min(150, n_cells), random_state=rng_seed)

    profile = [
        _mean_intensity_at_shift(img, sample, pixel_size_um, s, rng_seed) for s in PROFILE_SHIFTS_UM
    ]
    # shift=0 must be the maximum, with a small tolerance for sampling noise
    # (cells resampled at each shift distance, so exact monotonicity isn't
    # guaranteed even under perfect alignment).
    peak_at_zero = profile[0] >= max(profile[1:]) * 0.98

    return {
        "n_cells": n_cells,
        "n_cell_boundaries": n_cell_boundaries,
        "n_nucleus_boundaries": n_nucleus_boundaries,
        "coordinate_systems": list(sdata.coordinate_systems),
        "intensity_profile_by_shift_um": [round(v, 1) for v in profile],
        "peak_at_zero_shift": peak_at_zero,
    }


def build_coordinate_validation_report(
    spatialdata_root: Path,
    format_versions_path: Path,
    output_path: Path,
    project_root: Path,
) -> dict:
    import csv

    if not spatialdata_root.is_dir():
        raise PipelineError(
            f"'{spatialdata_root}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    with format_versions_path.open(newline="") as fh:
        pixel_sizes = {
            row["section_id"]: float(row["physical_size_x_um"])
            for row in csv.DictReader(fh, delimiter="\t")
        }

    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        raise PipelineError(f"No .zarr stores found under '{spatialdata_root}'.")

    if output_path.exists():
        output_path.unlink()
    writer = InventoryWriter(output_path, project_root=project_root, fields=VALIDATION_FIELDS)

    failures = []
    for zarr_path in zarr_paths:
        section_id = zarr_path.stem
        pixel_size = pixel_sizes.get(section_id)
        if pixel_size is None:
            raise PipelineError(
                f"No pixel size recorded for '{section_id}' in '{format_versions_path}'."
            )

        result = validate_section(zarr_path, pixel_size)
        passed = (
            result["peak_at_zero_shift"]
            and result["coordinate_systems"] == ["global"]
            and result["n_cells"] == result["n_cell_boundaries"]
        )
        status = "PASS" if passed else "FAIL"
        if not passed:
            failures.append(section_id)

        writer.write_row(
            section_id=section_id,
            n_cells=result["n_cells"],
            n_cell_boundaries=result["n_cell_boundaries"],
            n_nucleus_boundaries=result["n_nucleus_boundaries"],
            coordinate_systems=";".join(result["coordinate_systems"]),
            intensity_profile_by_shift_um=";".join(
                str(v) for v in result["intensity_profile_by_shift_um"]
            ),
            peak_at_zero_shift=result["peak_at_zero_shift"],
            status=status,
        )

    if failures:
        raise PipelineError(
            f"{len(failures)} section(s) failed coordinate validation: {failures}. See '{output_path}'."
        )

    return {
        "sections_validated": len(zarr_paths),
        "sections_passed": len(zarr_paths) - len(failures),
    }
