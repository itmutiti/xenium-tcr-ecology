"""Unit tests for xenium_tcr_ecology.clone_ecology.state_composition (`11_clone_spatial_descriptors/01_compute_clone_cell_state_composition.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.clone_ecology.state_composition import (
    STATE_ORDER,
    compute_shannon_entropy,
    compute_state_composition_with_ci,
)


class TestComputeStateCompositionWithCi:
    def test_all_one_state_gives_fraction_one_and_ci_below_one(self):
        states = pd.Series(["Cytotoxic"] * 5)
        result = compute_state_composition_with_ci(states)
        assert result["Cytotoxic"]["fraction"] == 1.0
        assert result["Cytotoxic"]["ci_low"] < 1.0
        assert result["Cytotoxic"]["ci_high"] == 1.0
        assert result["Exhausted"]["fraction"] == 0.0

    def test_mixed_states_fractions_sum_to_one(self):
        states = pd.Series(["Cytotoxic", "Cytotoxic", "Exhausted", "Cycling"])
        result = compute_state_composition_with_ci(states)
        total = sum(result[s]["fraction"] for s in STATE_ORDER)
        assert np.isclose(total, 1.0)
        assert result["Cytotoxic"]["fraction"] == 0.5

    def test_ci_widens_for_smaller_n_at_same_fraction(self):
        small = compute_state_composition_with_ci(pd.Series(["Cytotoxic", "Exhausted"]))
        large = compute_state_composition_with_ci(pd.Series(["Cytotoxic", "Exhausted"] * 20))
        small_width = small["Cytotoxic"]["ci_high"] - small["Cytotoxic"]["ci_low"]
        large_width = large["Cytotoxic"]["ci_high"] - large["Cytotoxic"]["ci_low"]
        assert small_width > large_width


class TestComputeShannonEntropy:
    def test_single_state_has_zero_entropy(self):
        assert compute_shannon_entropy([1.0, 0.0, 0.0]) == 0.0

    def test_uniform_two_states_has_entropy_one_bit(self):
        assert np.isclose(compute_shannon_entropy([0.5, 0.5]), 1.0)

    def test_uniform_four_states_has_entropy_two_bits(self):
        assert np.isclose(compute_shannon_entropy([0.25, 0.25, 0.25, 0.25]), 2.0)

    def test_zero_probability_states_contribute_nothing(self):
        result = compute_shannon_entropy([0.5, 0.5, 0.0, 0.0])
        assert np.isclose(result, 1.0)
