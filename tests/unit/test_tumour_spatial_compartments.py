"""Unit tests for xenium_tcr_ecology.tumour.spatial_compartments (`07_tumour_epithelium_characterisation/07_define_invasive_front_and_compartments.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.tumour.spatial_compartments import assign_compartment


class TestAssignCompartment:
    def test_deep_inside_is_tumour_core(self):
        result = assign_compartment(pd.Series([-10.0]), band_width_um=3.0)
        assert result.iloc[0] == "tumour_core"

    def test_just_inside_is_inner_margin(self):
        result = assign_compartment(pd.Series([-1.0]), band_width_um=3.0)
        assert result.iloc[0] == "inner_margin"

    def test_just_outside_is_outer_margin(self):
        result = assign_compartment(pd.Series([1.0]), band_width_um=3.0)
        assert result.iloc[0] == "outer_margin"

    def test_far_outside_is_distal_stroma(self):
        result = assign_compartment(pd.Series([100.0]), band_width_um=3.0)
        assert result.iloc[0] == "distal_stroma"

    def test_boundary_values_are_assigned_consistently(self):
        # Exactly at -band_width -> core (<=); exactly at 0 -> inner_margin
        # (<=0); exactly at +band_width -> outer_margin (<=).
        result = assign_compartment(pd.Series([-3.0, 0.0, 3.0]), band_width_um=3.0)
        assert list(result) == ["tumour_core", "inner_margin", "outer_margin"]

    def test_nan_distance_gives_nan_compartment(self):
        result = assign_compartment(pd.Series([np.nan, -10.0]), band_width_um=3.0)
        assert pd.isna(result.iloc[0])
        assert result.iloc[1] == "tumour_core"

    def test_narrower_band_reclassifies_shallow_interior_from_core_to_margin(self):
        # A cell at -2.0 is "core" under a 1.5um band but "inner_margin"
        # under a 3.0um band -- the whole point of the sensitivity check.
        wide = assign_compartment(pd.Series([-2.0]), band_width_um=3.0)
        narrow = assign_compartment(pd.Series([-2.0]), band_width_um=1.5)
        assert wide.iloc[0] == "inner_margin"
        assert narrow.iloc[0] == "tumour_core"
