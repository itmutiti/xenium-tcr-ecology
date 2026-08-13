"""Unit tests for xenium_tcr_ecology.io.xenium_reader (`03_spatialdata_import/01_import_each_section_to_spatialdata.py`).

_build_polygons is pure, fast logic and tested directly with synthetic
data. The full import_section() pipeline (h5/parquet/OME-TIFF parsing,
spatialdata model construction) is exercised against real project data in
test_real_project_import below, and was additionally validated manually
during development via a physical-alignment check (DAPI intensity ~10x
higher at real cell centroids than random background pixels) -- not
repeated here since it requires loading the full image into memory.
"""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.io.xenium_reader import _build_polygons


class TestBuildPolygons:
    def test_builds_one_polygon_per_cell(self):
        boundaries = pd.DataFrame(
            {
                "cell_id": ["c1", "c1", "c1", "c1", "c2", "c2", "c2"],
                "vertex_x": [0.0, 1.0, 1.0, 0.0, 5.0, 6.0, 5.5],
                "vertex_y": [0.0, 0.0, 1.0, 1.0, 5.0, 5.0, 6.0],
            }
        )
        gdf = _build_polygons(boundaries)
        assert set(gdf.index) == {"c1", "c2"}
        assert gdf.loc["c1", "geometry"].area == pytest.approx(1.0)

    def test_excludes_degenerate_boundaries(self):
        """A cell with fewer than 3 vertices cannot form a polygon -- must
        be excluded, not silently zero-filled or crashed on."""
        boundaries = pd.DataFrame(
            {
                "cell_id": ["c1", "c1", "c1", "c_bad", "c_bad"],
                "vertex_x": [0.0, 1.0, 0.5, 9.0, 9.0],
                "vertex_y": [0.0, 0.0, 1.0, 9.0, 9.0],
            }
        )
        gdf = _build_polygons(boundaries)
        assert "c_bad" not in gdf.index
        assert "c1" in gdf.index


class TestRealProjectImport:
    def test_real_section_imports_and_aligns(self):
        """Exercises the real reader against a real standardised section
        from this project's Raw Data Ingestion output, and repeats the physical
        alignment check (DAPI intensity at real cell centroids vs random
        background) as an automated regression test, not just a one-off
        manual check during development."""
        import numpy as np

        from xenium_tcr_ecology.infra.paths import find_project_root
        from xenium_tcr_ecology.io.xenium_reader import import_section

        project_root = find_project_root()
        section_dir = project_root / "data" / "standardised" / "P01_run1"
        if not section_dir.is_dir():
            pytest.skip("data/standardised/P01_run1 not present in this environment")

        sdata = import_section(section_dir, "P01_run1")

        assert sdata["table"].n_obs == sdata["cell_boundaries"].shape[0]
        assert sdata.coordinate_systems == ["global"]  # every element under one unified system

        img = sdata["morphology_mip"].compute().to_numpy()[0]
        cells = sdata["table"].obs.sample(20, random_state=42)
        px_size = 0.2125
        intensities = [
            img[int(row.y_centroid / px_size), int(row.x_centroid / px_size)]
            for _, row in cells.iterrows()
            if 0 <= int(row.y_centroid / px_size) < img.shape[0]
            and 0 <= int(row.x_centroid / px_size) < img.shape[1]
        ]
        rng = np.random.default_rng(42)
        bg = img[rng.integers(0, img.shape[0], 500), rng.integers(0, img.shape[1], 500)]

        assert np.mean(intensities) > 3 * np.mean(bg), (
            "Real cell centroids should sit on substantially higher DAPI signal than random "
            "background if the pixel<->micron coordinate alignment is correct."
        )
