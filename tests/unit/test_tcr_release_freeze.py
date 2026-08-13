"""Unit tests for xenium_tcr_ecology.tcr.release_freeze (`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.tcr.release_freeze import build_go_no_go_decision


class TestBuildGoNoGoDecision:
    def test_status_is_conditional_go(self):
        clone_metadata = pd.DataFrame({"n_singlet_cells": [5, 0, 3]})
        ascertainment = pd.DataFrame({"intended_patient_identified": [True, False, True]})
        result = build_go_no_go_decision(clone_metadata, ascertainment)
        assert result["status"] == "CONDITIONAL GO"

    def test_caveats_are_non_empty_and_reference_real_figures(self):
        clone_metadata = pd.DataFrame({"n_singlet_cells": [5, 0, 3, 0]})
        ascertainment = pd.DataFrame({"intended_patient_identified": [True, False, True, True]})
        result = build_go_no_go_decision(clone_metadata, ascertainment)
        assert len(result["caveats"]) >= 3
        # 1/4 not identified -> 25.0%
        assert any("25.0%" in c for c in result["caveats"])
        # 2/4 clones excluded -> 50.0%
        assert any("50.0%" in c for c in result["caveats"])
