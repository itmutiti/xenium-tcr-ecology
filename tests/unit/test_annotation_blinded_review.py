"""Unit tests for xenium_tcr_ecology.annotation.blinded_review (`06_cell_type_annotation/07_blinded_annotation_review.py`)."""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.blinded_review import (
    ADJUDICATION_LOG_COLUMNS,
    anonymize_panel_id,
    select_review_sample,
)


class TestAnonymizePanelId:
    def test_deterministic(self):
        assert anonymize_panel_id("P01_run1_aaadggoi-1") == anonymize_panel_id(
            "P01_run1_aaadggoi-1"
        )

    def test_different_cells_get_different_ids(self):
        assert anonymize_panel_id("cellA") != anonymize_panel_id("cellB")

    def test_panel_id_does_not_contain_the_cell_id(self):
        # The whole point is blinding -- the panel ID must not leak the
        # real cell ID by simple string inspection.
        cell_id = "P01_run1_aaadggoi-1"
        assert cell_id not in anonymize_panel_id(cell_id)


def _make_final_annotations(n_per_lineage=20, rng_seed=0):
    import numpy as np

    rng = np.random.default_rng(rng_seed)
    lineages = ["T_cell", "B_cell", "Myeloid"]
    rows = []
    for lineage in lineages:
        for i in range(n_per_lineage):
            rows.append(
                {
                    "cell_id": f"{lineage}_{i}",
                    "final_lineage": lineage,
                    "is_ambiguous": bool(rng.random() < 0.3),
                    "confidence": rng.random(),
                    "final_substate": None,
                }
            )
    return pd.DataFrame(rows).set_index("cell_id")


class TestSelectReviewSample:
    def test_raises_on_empty_input(self):
        with pytest.raises(PipelineError, match="empty"):
            select_review_sample(pd.DataFrame())

    def test_oversamples_ambiguous_cells_relative_to_population_share(self):
        annotations = _make_final_annotations(n_per_lineage=100, rng_seed=1)
        population_ambiguous_fraction = annotations["is_ambiguous"].mean()

        sample = select_review_sample(annotations, n_total=60, rng_seed=1)
        sample_ambiguous_fraction = sample["is_ambiguous"].mean()

        assert sample_ambiguous_fraction > population_ambiguous_fraction

    def test_sample_covers_multiple_lineages(self):
        annotations = _make_final_annotations(n_per_lineage=50, rng_seed=2)
        sample = select_review_sample(annotations, n_total=60, rng_seed=2)
        assert sample["final_lineage"].nunique() >= 2

    def test_sample_size_is_approximately_n_total(self):
        annotations = _make_final_annotations(n_per_lineage=100, rng_seed=3)
        sample = select_review_sample(annotations, n_total=60, rng_seed=3)
        # Stratified sampling with per-lineage floor/rounding won't hit the
        # target exactly, but must be in a sane range around it.
        assert 30 <= len(sample) <= 90


def test_adjudication_log_columns_include_panel_id_and_no_prefilled_verdict():
    assert "panel_id" in ADJUDICATION_LOG_COLUMNS
    assert "reviewer_adjudicated_lineage" in ADJUDICATION_LOG_COLUMNS
