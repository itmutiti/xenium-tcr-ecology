"""Unit tests for xenium_tcr_ecology.release.redesigned_main_figures (`17_statistical_closure_and_release/11_build_redesigned_manuscript_figures.py`)."""

from __future__ import annotations

import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.release.redesigned_main_figures import (
    build_figure_1_framework_generalisation,
    build_figure_barrier_topology_with_ablation,
    build_figure_hpv_consolidated,
    build_figure_variance_partition_with_sensitivity,
)

NULL_MODEL_COLS = [
    "pvalue_constrained_permutation",
    "pvalue_degree_preserving",
    "pvalue_graph_preserving",
]


def _write_calibration_results(path, n_replicates: int = 5) -> None:
    rows = []
    for replicate in range(n_replicates):
        for effect_size in [0.0, 0.1, 0.5, 1.0]:
            rows.append(
                {
                    "replicate": replicate,
                    "effect_size": effect_size,
                    **{c: 0.5 - effect_size * 0.4 for c in NULL_MODEL_COLS},
                }
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(path)


class TestBuildFigure1FrameworkGeneralisation:
    def test_real_composite_figure_written(self, tmp_path):
        _write_calibration_results(
            tmp_path / "reports" / "graphs" / "null_model_calibration.parquet"
        )
        _write_calibration_results(
            tmp_path / "data" / "derived" / "framework_generalisation_results.parquet"
        )
        _write_calibration_results(
            tmp_path
            / "data"
            / "derived"
            / "framework_generalisation_results_second_dataset.parquet"
        )

        result = build_figure_1_framework_generalisation(tmp_path)
        output_path = (
            tmp_path
            / "reports"
            / "manuscript_figures"
            / "framework_generalisation_three_tumour_types.pdf"
        )
        assert output_path.is_file()
        assert output_path.stat().st_size > 0
        assert result["output_path"] == str(output_path)

    def test_real_missing_source_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_figure_1_framework_generalisation(tmp_path)


class TestBuildFigureVariancePartitionWithSensitivity:
    def _write_partition(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "component": ["patient", "identity", "context"],
                "variance": [0.3, 0.2, 0.5],
                "proportion": [0.3, 0.2, 0.5],
                "ci_low": [0.05, 0.05, 0.3],
                "ci_high": [0.5, 0.4, 0.7],
            }
        ).to_parquet(path)

    def test_real_composite_figure_written(self, tmp_path):
        self._write_partition(tmp_path / "data" / "derived" / "variance_partition_results.parquet")
        self._write_partition(
            tmp_path
            / "data"
            / "derived"
            / "variance_partition_sensitivity_excluding_cycling.parquet"
        )

        result = build_figure_variance_partition_with_sensitivity(tmp_path)
        output_path = (
            tmp_path / "reports" / "manuscript_figures" / "variance_partition_with_sensitivity.pdf"
        )
        assert output_path.is_file()
        assert result["output_path"] == str(output_path)

    def test_real_missing_source_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_figure_variance_partition_with_sensitivity(tmp_path)


class TestBuildFigureBarrierTopologyWithAblation:
    def test_real_composite_figure_written(self, tmp_path):
        model_path = tmp_path / "data" / "derived" / "barrier_topology_model_results.parquet"
        model_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "covariate": [
                    "fibroblast_barrier_fraction",
                    "suppressive_myeloid_barrier_fraction",
                ],
                "estimate": [0.05, -0.34],
                "se": [0.07, 0.11],
                "ci_low": [-0.08, -0.56],
                "ci_high": [0.18, -0.13],
            }
        ).to_parquet(model_path)

        ablation_path = tmp_path / "data" / "derived" / "barrier_covariate_ablation.parquet"
        pd.DataFrame(
            {
                "step": [
                    "barrier_only",
                    "state_block_only",
                    "niche_block_only",
                    "+ niche_archetype_4_fraction",
                    "full_state_and_niche",
                ],
                "n_adjustment_covariates": [0, 4, 5, 1, 9],
                "estimate": [0.06, -0.05, -0.16, -0.21, -0.34],
                "se": [0.13, 0.14, 0.11, 0.12, 0.11],
                "ci_low": [-0.2, -0.32, -0.38, -0.44, -0.56],
                "ci_high": [0.32, 0.22, 0.05, 0.02, -0.13],
                "p_value": [0.66, 0.72, 0.13, 0.08, 0.0015],
                "marginal_r2": [0.001, 0.02, 0.05, 0.09, 0.618],
            }
        ).to_parquet(ablation_path)

        literature_path = tmp_path / "data" / "derived" / "literature_benchmark_results.parquet"
        pd.DataFrame(
            {
                "covariate": [
                    "fibroblast_barrier_fraction",
                    "suppressive_myeloid_barrier_fraction",
                ],
                "project_raw_r": [0.019, -0.079],
                "project_raw_pvalue": [0.8, 0.33],
                "published_r": [None, -0.48],
                "published_pvalue": [None, 1e-10],
                "published_citation": [None, "Grout et al. 2022"],
                "direction_concordant": [None, True],
                "magnitude_ratio": [None, 0.165],
            }
        ).to_parquet(literature_path)

        result = build_figure_barrier_topology_with_ablation(tmp_path)
        output_path = (
            tmp_path / "reports" / "manuscript_figures" / "barrier_topology_with_ablation.pdf"
        )
        assert output_path.is_file()
        assert result["output_path"] == str(output_path)

    def test_real_missing_source_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_figure_barrier_topology_with_ablation(tmp_path)


class TestBuildFigureHpvConsolidated:
    def test_real_composite_figure_written(self, tmp_path):
        composition_path = (
            tmp_path / "data" / "derived" / "hpv_composition_comparison_results.parquet"
        )
        composition_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "lineage": ["T_cell", "B_cell"],
                "n_positive": [4, 4],
                "n_negative": [4, 4],
                "median_positive": [0.1, 0.03],
                "median_negative": [0.08, 0.03],
                "pvalue": [0.34, 0.68],
                "pvalue_bh": [0.45, 0.88],
            }
        ).to_parquet(composition_path)

        structure_path = tmp_path / "data" / "derived" / "hpv_structure_comparison_results.parquet"
        pd.DataFrame(
            {
                "outcome_domain": ["ecosystem", "clone_structure"],
                "category": ["Tumour niche", "ecological_structure_score"],
                "metric": ["abundance", "ecological_structure_score"],
                "n_positive": [4, 4],
                "n_negative": [4, 4],
                "median_positive": [0.2, 0.1],
                "median_negative": [0.18, 0.09],
                "pvalue": [0.49, 0.6],
                "pvalue_bh": [0.89, 0.9],
            }
        ).to_parquet(structure_path)

        hpv_claims_path = tmp_path / "tables" / "Table_10_hpv_claim_strength.tsv"
        hpv_claims_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "claim_id": ["hpv_discordance_P01", "hpv_probe_positive_untested_P17"],
                "claim": [
                    "Patient P01's clinical p16 status (Positive) is discordant with its real Xenium HPV16 E6/E7 probe signal (0.0025 of cells).",
                    "Patient P17 shows real detectable HPV16 E6/E7 signal (0.1948 of cells) despite no clinical p16 test.",
                ],
            }
        ).to_csv(hpv_claims_path, sep="\t", index=False)

        hpv_metadata_path = tmp_path / "metadata" / "hpv_status_validated.tsv"
        hpv_metadata_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "patient_id": ["P01", "P03", "P09", "P17"],
                "clinical_p16_status": ["Positive", "Positive", "Positive", "Not Tested"],
                "has_hpv_probe_coverage": [True, True, False, True],
                "hpv_e6_e7_probe_positive_fraction": [0.0025, 0.42, None, 0.1948],
                "validated_hpv_status": [
                    "discordant_clinical_positive_probe_negative",
                    "confirmed_positive",
                    "presumed_negative_unverifiable",
                    "probe_positive_clinically_untested",
                ],
            }
        ).to_csv(hpv_metadata_path, sep="\t", index=False)

        result = build_figure_hpv_consolidated(tmp_path)
        output_path = (
            tmp_path
            / "reports"
            / "manuscript_figures"
            / "hpv_composition_structure_and_discordance_qc.pdf"
        )
        assert output_path.is_file()
        assert result["output_path"] == str(output_path)

    def test_real_missing_source_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_figure_hpv_consolidated(tmp_path)
