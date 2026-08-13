"""Unit tests for xenium_tcr_ecology.qc.qc_release_report (`04_quality_control/09_generate_qc_release_report.py`)."""

from __future__ import annotations

import pandas as pd
import pytest
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.qc.qc_release_report import (
    build_go_no_go_decision,
    build_qc_release_report,
)


def _write_fixture_inputs(project_root):
    (project_root / "data" / "derived").mkdir(parents=True, exist_ok=True)
    (project_root / "reports" / "qc" / "spatial_artifact_masks").mkdir(parents=True, exist_ok=True)
    (project_root / "config").mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        {
            "section_id": ["S1_run1", "S1_run1", "S2_run1"],
            "transcript_counts": [100, 150, 200],
            "n_genes_detected": [20, 25, 30],
        }
    ).to_parquet(project_root / "data" / "derived" / "cell_qc_metrics.parquet")

    pd.DataFrame(
        {
            "section_id": ["S1_run1", "S2_run1"],
            "fraction_qv_below_20": [0.1, 0.12],
            "fraction_overlaps_nucleus": [0.4, 0.42],
        }
    ).to_parquet(project_root / "data" / "derived" / "transcript_qc_metrics.parquet")

    pd.DataFrame(
        {
            "section_id": ["S1_run1", "S1_run1", "S2_run1"],
            "fov_name": ["A1", "A2", "B1"],
            "flagged_artifact_candidate": [False, True, False],
        }
    ).to_csv(
        project_root / "reports" / "qc" / "spatial_artifact_masks" / "fov_artifact_candidates.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        {
            "section_id": ["S1_run1", "S2_run1"],
            "fraction_invalid_cell_polygon": [0.0, 0.01],
            "fraction_invalid_nucleus_polygon": [0.02, 0.03],
            "fraction_nucleus_not_contained": [0.0, 0.0],
        }
    ).to_parquet(project_root / "reports" / "qc" / "segmentation_review.parquet")

    with open(project_root / "config" / "qc_thresholds.yaml", "w") as f:
        yaml.dump({"active_profile": "standard"}, f)

    pd.DataFrame(
        {
            "section_id": ["S1_run1", "S1_run1", "S2_run1"],
            "qc_pass": [True, False, True],
            "excluded": [False, True, False],
        }
    ).to_csv(project_root / "data" / "derived" / "exclusion_log.tsv", sep="\t", index=False)

    pd.DataFrame(
        {
            "patient_id": ["S1"],
            "pseudobulk_pearson_r": [0.98],
            "flagged_discordant": [False],
        }
    ).to_csv(project_root / "data" / "derived" / "replicate_concordance.tsv", sep="\t", index=False)


class TestBuildGoNoGoDecision:
    def test_status_is_conditional_go(self, tmp_path):
        path = tmp_path / "replicate_concordance.tsv"
        pd.DataFrame({"patient_id": ["P1", "P2"], "flagged_discordant": [False, False]}).to_csv(
            path, sep="\t", index=False
        )
        decision = build_go_no_go_decision(path)
        assert decision["status"] == "CONDITIONAL GO"
        assert "4.04" in " ".join(decision["outstanding_workstreams"])
        assert "4.05" in " ".join(decision["outstanding_workstreams"])

    def test_counts_flagged_discordant_pairs(self, tmp_path):
        path = tmp_path / "replicate_concordance.tsv"
        pd.DataFrame(
            {"patient_id": ["P1", "P2", "P3"], "flagged_discordant": [False, True, False]}
        ).to_csv(path, sep="\t", index=False)
        decision = build_go_no_go_decision(path)
        assert decision["replicate_concordance_pairs_flagged"] == 1
        assert decision["replicate_concordance_flagged_patients"] == ["P2"]


class TestBuildQcReleaseReport:
    def test_writes_report_and_summary(self, tmp_path):
        _write_fixture_inputs(tmp_path)
        output_path = tmp_path / "reports" / "qc" / "QC_release_report.html"

        summary = build_qc_release_report(tmp_path, output_path)

        assert summary["status"] == "CONDITIONAL GO"
        assert summary["sections_processed"] == 2
        assert summary["n_cells_raw_total"] == 3
        assert output_path.is_file()
        assert (tmp_path / "reports" / "qc" / "qc_release_section_summary.parquet").is_file()
        assert "CONDITIONAL GO" in output_path.read_text()

    def test_raises_on_missing_input(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            build_qc_release_report(tmp_path, tmp_path / "out.html")
