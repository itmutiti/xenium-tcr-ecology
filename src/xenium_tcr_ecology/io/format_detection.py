"""Xenium output format detection and standard-reader compatibility check
(`03_spatialdata_import/00_detect_xenium_format_version.py`).

Verified against this dataset (not from
generic documentation): GEO's supplementary release for GSE300147 does NOT
include ``experiment.xenium`` or ``cells.zarr.zip``, both of which
spatialdata_io.xenium() requires unconditionally (it opens
``path / "experiment.xenium"`` with no fallback). The standard reader
therefore cannot be used on this dataset -- a custom reader (`03_spatialdata_import/01_import_each_section_to_spatialdata.py`)
is required, built directly from the 6 files GEO does publish:
cell_feature_matrix.h5, cells.parquet.gz, cell_boundaries.parquet.gz,
nucleus_boundaries.parquet.gz, transcripts.parquet.gz, morphology.ome.tif.gz.

Also extracts the authoritative pixel-to-micron scale from each section's
OME-XML metadata (PhysicalSizeX/Y/Z), reading directly from the gzip stream
without a full decompression, and confirms it is consistent across all
sections -- this is the coordinate transform `03_spatialdata_import/01_import_each_section_to_spatialdata.py` needs to align the
morphology image (pixel space) with cell/transcript coordinates (already in
micron space, confirmed separately against cell centroid ranges).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import tifffile

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter

STANDARD_READER_REQUIRED_FILES = ["experiment.xenium", "cells.zarr.zip"]
CANONICAL_FILES = [
    "cell_boundaries.parquet.gz",
    "cell_feature_matrix.h5",
    "cells.parquet.gz",
    "morphology.ome.tif.gz",
    "nucleus_boundaries.parquet.gz",
    "transcripts.parquet.gz",
]

FORMAT_FIELDS = [
    "section_id",
    "standard_reader_compatible",
    "missing_for_standard_reader",
    "size_x_px",
    "size_y_px",
    "size_z",
    "physical_size_x_um",
    "physical_size_y_um",
    "physical_size_z_um",
    "channel_name",
]


def extract_ome_pixel_metadata(morphology_gz_path: Path) -> dict:
    with gzip.open(morphology_gz_path, "rb") as f:
        tf = tifffile.TiffFile(f)
        shape = tf.series[0].shape  # (Z, Y, X) for this dataset, verified
        ome_xml = tf.ome_metadata

    import xml.etree.ElementTree as ET

    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    root = ET.fromstring(ome_xml)
    pixels = root.find(".//ome:Pixels", ns)
    channel = root.find(".//ome:Channel", ns)

    return {
        "size_z_from_series": shape[0],
        "size_y_from_series": shape[1],
        "size_x_from_series": shape[2],
        "size_x_px": int(pixels.get("SizeX")),
        "size_y_px": int(pixels.get("SizeY")),
        "size_z": int(pixels.get("SizeZ")),
        "physical_size_x_um": float(pixels.get("PhysicalSizeX")),
        "physical_size_y_um": float(pixels.get("PhysicalSizeY")),
        "physical_size_z_um": float(pixels.get("PhysicalSizeZ")),
        "channel_name": channel.get("Name") if channel is not None else None,
    }


def check_standard_reader_compatibility(section_dir: Path) -> dict:
    present = {f.name for f in section_dir.iterdir()}
    missing = [f for f in STANDARD_READER_REQUIRED_FILES if f not in present]
    return {"compatible": len(missing) == 0, "missing": missing}


def build_format_report(standardised_root: Path, output_path: Path, project_root: Path) -> dict:
    if not standardised_root.is_dir():
        raise PipelineError(
            f"Standardised data directory not found: '{standardised_root}'. Run `02_raw_data_ingestion/05_standardise_sample_directory_layout.py` first."
        )

    section_dirs = sorted(p for p in standardised_root.iterdir() if p.is_dir())
    if not section_dirs:
        raise PipelineError(f"No section directories found under '{standardised_root}'.")

    if output_path.exists():
        output_path.unlink()
    writer = InventoryWriter(output_path, project_root=project_root, fields=FORMAT_FIELDS)

    pixel_sizes_seen = set()
    incompatible_count = 0
    for section_dir in section_dirs:
        section_id = section_dir.name
        for canonical in CANONICAL_FILES:
            if not (section_dir / canonical).exists():
                raise PipelineError(
                    f"Section '{section_id}' is missing canonical file '{canonical}' -- `02_raw_data_ingestion/05_standardise_sample_directory_layout.py` should have prevented this."
                )

        compat = check_standard_reader_compatibility(section_dir)
        if not compat["compatible"]:
            incompatible_count += 1

        ome = extract_ome_pixel_metadata(section_dir / "morphology.ome.tif.gz")
        pixel_sizes_seen.add(
            (ome["physical_size_x_um"], ome["physical_size_y_um"], ome["physical_size_z_um"])
        )

        writer.write_row(
            section_id=section_id,
            standard_reader_compatible=compat["compatible"],
            missing_for_standard_reader=";".join(compat["missing"]),
            size_x_px=ome["size_x_px"],
            size_y_px=ome["size_y_px"],
            size_z=ome["size_z"],
            physical_size_x_um=ome["physical_size_x_um"],
            physical_size_y_um=ome["physical_size_y_um"],
            physical_size_z_um=ome["physical_size_z_um"],
            channel_name=ome["channel_name"],
        )

    return {
        "sections_checked": len(section_dirs),
        "standard_reader_incompatible": incompatible_count,
        "distinct_pixel_size_um_combos": len(pixel_sizes_seen),
        "pixel_sizes_um": sorted(pixel_sizes_seen),
    }
