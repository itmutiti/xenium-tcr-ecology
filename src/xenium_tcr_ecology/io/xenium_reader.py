"""Custom Xenium section reader (`03_spatialdata_import/01_import_each_section_to_spatialdata.py`).

spatialdata_io.xenium() cannot be used on this dataset (`03_spatialdata_import/00_detect_xenium_format_version.py`: GEO
publishes neither experiment.xenium nor cells.zarr.zip). This module builds
a SpatialData object directly from the 6 files GEO does publish, using
spatialdata's element models (Image2DModel, ShapesModel, PointsModel,
TableModel) rather than any private/internal spatialdata_io machinery.

Coordinate systems, verified against the data (`03_spatialdata_import/00_detect_xenium_format_version.py` and
manual inspection):
  - cells.parquet's x_centroid/y_centroid, and transcripts.parquet's
    x_location/y_location, are already in microns.
  - cell_boundaries.parquet / nucleus_boundaries.parquet vertex_x/vertex_y
    are in the same micron space (same coordinate frame as the above).
  - morphology.ome.tif is in pixel space; OME-XML PhysicalSizeX/Y (0.2125
    um/pixel for every section in this cohort, confirmed in `03_spatialdata_import/00_detect_xenium_format_version.py`) is
    the scale factor to align it with the micron space above.

The morphology image is a focus stack (11 Z-planes at 3.0 um spacing) of a
single 2D tissue section, not a true 3D volume of distinct cell positions
(cells/transcripts have one x,y per cell; only transcripts carry a z, which
is a per-molecule z-offset within the section, not evidence of 3D tissue
structure). A maximum-intensity projection (MIP) across Z is used for the
2D image element -- the same convention 10x's own "morphology_mip" output
represents, which this GEO release does not separately publish. This is a
documented design decision, not an incidental default.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import tifffile
from geopandas import GeoDataFrame
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Identity, Scale, set_transformation

from xenium_tcr_ecology.infra.exceptions import PipelineError

# "global" (not a custom name like "microns"): spatialdata's element models
# register an implicit Identity transform under "global" by default at
# parse() time. Using any other coordinate-system name for shapes/points
# leaves them double-registered (their own default "global" entry, plus
# the custom one) while an image given an explicit `transformations` dict
# at parse time is registered only under the name given -- verified
# directly against actual output: naming this "microns" left the image
# absent from the "global" system entirely, so anything querying elements
# in "global" would silently miss it. Standardising on "global" for every
# element's transformation collapses this to one single, consistent system.
COORDINATE_SYSTEM = "global"


def _read_gzipped_parquet(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rb") as f:
        return pd.read_parquet(io.BytesIO(f.read()))


def _build_polygons(boundaries: pd.DataFrame) -> GeoDataFrame:
    """boundaries: long-format DataFrame with columns cell_id, vertex_x, vertex_y,
    one row per polygon vertex, in order. Builds one Polygon per cell_id."""
    polygons = {}
    for cell_id, group in boundaries.groupby("cell_id", sort=False):
        coords = list(zip(group["vertex_x"].to_numpy(), group["vertex_y"].to_numpy()))
        if len(coords) < 3:
            continue  # degenerate boundary, cannot form a polygon -- excluded, not silently zero-filled
        polygons[cell_id] = Polygon(coords)
    gdf = GeoDataFrame({"geometry": list(polygons.values())}, index=list(polygons.keys()))
    return gdf


def _read_morphology_mip(
    morphology_gz_path: Path, physical_size_x_um: float, physical_size_y_um: float
):
    with gzip.open(morphology_gz_path, "rb") as f:
        tf = tifffile.TiffFile(f)
        stack = tf.series[
            0
        ].asarray()  # (Z, Y, X), verified in `03_spatialdata_import/00_detect_xenium_format_version.py`

    mip = stack.max(axis=0)  # (Y, X) -- see module docstring for why MIP, not a single plane
    image_data = mip[
        np.newaxis, :, :
    ]  # (C=1, Y, X), the shape Image2DModel expects with dims "cyx"

    image = Image2DModel.parse(
        image_data,
        dims=("c", "y", "x"),
        transformations={
            COORDINATE_SYSTEM: Scale([physical_size_x_um, physical_size_y_um], axes=("x", "y"))
        },
    )
    return image


def import_section(section_dir: Path, section_id: str) -> SpatialData:
    required = [
        "cell_feature_matrix.h5",
        "cells.parquet.gz",
        "cell_boundaries.parquet.gz",
        "nucleus_boundaries.parquet.gz",
        "transcripts.parquet.gz",
        "morphology.ome.tif.gz",
    ]
    missing = [f for f in required if not (section_dir / f).is_file()]
    if missing:
        raise PipelineError(f"Section '{section_id}' at '{section_dir}' is missing: {missing}")

    # --- expression table ---
    adata = sc.read_10x_h5(section_dir / "cell_feature_matrix.h5")
    adata.var_names_make_unique()

    # --- cell metadata, joined onto the table ---
    cells = _read_gzipped_parquet(section_dir / "cells.parquet.gz").set_index("cell_id")
    missing_cells = set(adata.obs_names) - set(cells.index)
    if missing_cells:
        raise PipelineError(
            f"Section '{section_id}': {len(missing_cells)} cell(s) in the expression matrix have no "
            f"entry in cells.parquet.gz -- e.g. {sorted(missing_cells)[:5]}"
        )
    adata.obs = adata.obs.join(cells, how="left")
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    adata.obs["region"] = "cell_boundaries"
    adata.obs["cell_id"] = adata.obs_names

    table = TableModel.parse(
        adata, region="cell_boundaries", region_key="region", instance_key="cell_id"
    )

    # --- shapes: cell and nucleus boundaries ---
    cell_boundaries_raw = _read_gzipped_parquet(section_dir / "cell_boundaries.parquet.gz")
    nucleus_boundaries_raw = _read_gzipped_parquet(section_dir / "nucleus_boundaries.parquet.gz")
    cell_shapes = ShapesModel.parse(_build_polygons(cell_boundaries_raw))
    nucleus_shapes = ShapesModel.parse(_build_polygons(nucleus_boundaries_raw))
    set_transformation(cell_shapes, Identity(), COORDINATE_SYSTEM)
    set_transformation(nucleus_shapes, Identity(), COORDINATE_SYSTEM)

    # --- points: transcripts ---
    transcripts = _read_gzipped_parquet(section_dir / "transcripts.parquet.gz")
    points = PointsModel.parse(
        transcripts, coordinates={"x": "x_location", "y": "y_location", "z": "z_location"}
    )
    set_transformation(points, Identity(), COORDINATE_SYSTEM)

    # --- image: morphology MIP, scaled from pixel to micron space ---
    import xml.etree.ElementTree as ET

    with gzip.open(section_dir / "morphology.ome.tif.gz", "rb") as f:
        ome_xml = tifffile.TiffFile(f).ome_metadata
    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    pixels_el = ET.fromstring(ome_xml).find(".//ome:Pixels", ns)
    px_x, px_y = float(pixels_el.get("PhysicalSizeX")), float(pixels_el.get("PhysicalSizeY"))
    image = _read_morphology_mip(section_dir / "morphology.ome.tif.gz", px_x, px_y)

    sdata = SpatialData(
        images={"morphology_mip": image},
        shapes={"cell_boundaries": cell_shapes, "nucleus_boundaries": nucleus_shapes},
        points={"transcripts": points},
        tables={"table": table},
    )
    return sdata
