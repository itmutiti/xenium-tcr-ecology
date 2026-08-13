"""Unit tests for xenium_tcr_ecology.clone_ecology.clone_atlas (`13_clone_ecology_confirmatory_models/06_generate_clone_atlas.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_tcr_ecology.clone_ecology.clone_atlas import (
    build_clone_card_html,
    render_clone_thumbnail,
)


class TestRenderCloneThumbnail:
    def test_returns_a_real_embeddable_png_data_uri(self):
        rng = np.random.default_rng(0)
        section_x = rng.normal(size=100)
        section_y = rng.normal(size=100)
        clone_x = section_x[:5]
        clone_y = section_y[:5]
        result = render_clone_thumbnail(clone_x, clone_y, section_x, section_y)
        assert result.startswith("data:image/png;base64,")
        assert len(result) > 100  # a real, non-trivial image was encoded

    def test_handles_a_single_cell_clone_without_crashing(self):
        result = render_clone_thumbnail(
            np.array([1.0]), np.array([2.0]), np.array([1.0, 2.0, 3.0]), np.array([2.0, 3.0, 4.0])
        )
        assert result.startswith("data:image/png;base64,")


class TestBuildCloneCardHtml:
    def test_real_html_contains_clone_and_patient_identifiers(self):
        row = pd.Series(
            {
                "clone_id": "230322_TESTCLONE_TRA",
                "patient_id": "P01",
                "section_id": "P01_run1",
                "n_cells": 5,
                "ecological_structure_score": 0.5,
                "dominant_state": "Exhausted",
            }
        )
        result = build_clone_card_html(row, "data:image/png;base64,ABC123")
        assert "230322_TESTCLONE_TRA" in result
        assert "P01" in result
        assert "P01_run1" in result
        assert "ABC123" in result

    def test_html_special_characters_in_identifiers_are_escaped(self):
        row = pd.Series(
            {
                "clone_id": "clone<script>alert(1)</script>",
                "patient_id": "P01",
                "section_id": "P01_run1",
                "n_cells": 5,
                "ecological_structure_score": 0.5,
                "dominant_state": "Exhausted",
            }
        )
        result = build_clone_card_html(row, "data:image/png;base64,ABC123")
        assert "<script>" not in result
        assert "&lt;script&gt;" in result
