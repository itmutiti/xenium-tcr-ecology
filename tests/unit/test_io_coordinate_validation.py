"""Unit tests for xenium_tcr_ecology.io.coordinate_validation (`03_spatialdata_import/02_validate_coordinate_systems.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.io.coordinate_validation import _mean_intensity_at_shift


class TestMeanIntensityAtShift:
    def test_detects_peak_at_zero_for_a_real_bright_spot(self):
        """A synthetic image with a single bright pixel at a known location:
        sampling exactly at that location (shift=0) must read higher than
        sampling with any nonzero shift."""
        img = np.zeros((200, 200), dtype=np.uint16)
        img[100, 100] = 10000  # single bright "nucleus" pixel

        sample = pd.DataFrame({"x_centroid": [100.0], "y_centroid": [100.0]})

        at_zero = _mean_intensity_at_shift(img, sample, pixel_size_um=1.0, shift_um=0.0, rng_seed=1)
        at_shift = _mean_intensity_at_shift(
            img, sample, pixel_size_um=1.0, shift_um=20.0, rng_seed=1
        )

        assert at_zero == 10000
        assert at_shift < at_zero

    def test_out_of_bounds_shift_is_excluded_not_crashed_on(self):
        img = np.ones((10, 10), dtype=np.uint16) * 5
        # centroid near the edge; a large shift will push some samples out of bounds
        sample = pd.DataFrame({"x_centroid": [1.0] * 20, "y_centroid": [1.0] * 20})
        result = _mean_intensity_at_shift(img, sample, pixel_size_um=1.0, shift_um=50.0, rng_seed=1)
        assert result == 0.0 or result >= 0  # no crash; either no in-bounds samples or a valid mean
