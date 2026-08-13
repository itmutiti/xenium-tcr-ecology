"""Unit tests for xenium_tcr_ecology.external_checkpoint.freeze_decision (`12_external_checkpoint_validation/03_decide_freeze_or_revise.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.external_checkpoint.freeze_decision import (
    decide_freeze_or_revise,
    identify_flagged_states,
)


class TestDecideFreezeOrRevise:
    def test_both_conditions_met_gives_freeze(self):
        assert (
            decide_freeze_or_revise(all_programs_transfer=True, overall_spearman_rho=0.5)
            == "freeze"
        )

    def test_programs_not_transferring_gives_revise(self):
        assert (
            decide_freeze_or_revise(all_programs_transfer=False, overall_spearman_rho=0.9)
            == "revise"
        )

    def test_negative_rho_gives_revise(self):
        assert (
            decide_freeze_or_revise(all_programs_transfer=True, overall_spearman_rho=-0.2)
            == "revise"
        )

    def test_zero_rho_gives_revise(self):
        # rho == 0 is not a real positive direction -- must not freeze on a null signal.
        assert (
            decide_freeze_or_revise(all_programs_transfer=True, overall_spearman_rho=0.0)
            == "revise"
        )


class TestIdentifyFlaggedStates:
    def test_flags_states_above_threshold(self):
        df = pd.DataFrame({"state": ["A", "B", "C"], "rank_shift": [4, -3, 1]})
        result = identify_flagged_states(df, threshold=3)
        assert result == ["A", "B"]

    def test_no_states_flagged_when_all_small(self):
        df = pd.DataFrame({"state": ["A", "B"], "rank_shift": [1, -2]})
        result = identify_flagged_states(df, threshold=3)
        assert result == []
