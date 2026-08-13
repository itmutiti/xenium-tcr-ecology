"""Links every registered claim (`governance/analysis_registry.tsv`) to
its discovery, sensitivity, replicate and external-validation evidence,
with an overall grade (`16_external_validation_and_generalisation/07_generate_evidence_matrix.py`) -- the final synthesis milestone of External Validation and Generalisation.

**Complete coverage of the registered claim set:** each of the 11
analyses in `governance/analysis_registry.tsv` gets exactly one evidence
row here (checked structurally by
`test_real_every_registered_analysis_is_covered`), not just the
confirmatory Q1-Q3/HPV claims -- methods-validation gates
(`null_model_calibration_suite`, `tcr_false_positive_rate`,
`replicate_concordance`) and robustness/consistency checks
(`segmentation_robustness_check`, `external_checkpoint_directional_
consistency`) are claims too.

**Grading is not uniformly optimistic:** `overall_evidence_grade` ranges
from `strong` (multiple independent lines of evidence converge, e.g.
`q1_framework_generalisation`, `tcr_false_positive_rate`) through
`moderate` (significant, with a reported caveat, e.g.
`q3_barrier_topology_confirmatory`'s magnitude weaker than the closest
published comparator) down to `weak_exploratory` (`hpv_primary_contrast`
-- 0/25 BH-significant tests, only 2 exploratory leads, consistent with
`15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`'s
power simulation) and `supporting` (analyses that are themselves
evidence for another claim, not independent claims needing their own
full evidentiary stack, e.g. `q3_literature_benchmark`,
`segmentation_robustness_check`).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

VALID_GRADES = {"strong", "moderate", "weak_exploratory", "supporting"}

CLAIM_EVIDENCE_MATRIX: list[dict] = [
    {
        "claim_id": "q1_framework_generalisation",
        "discovery_evidence": "`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`: null-model calibration suite on `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s synthetic ground truth (10 synthetic patients).",
        "sensitivity_evidence": "n/a -- this claim IS itself the sensitivity/calibration check for all downstream spatial-association claims.",
        "replicate_evidence": "10 synthetic-patient replicates (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`) + 10 independent subsample replicates on a different dataset (`16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py`).",
        "external_validation_evidence": "`16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py`: 3/3 null models CI-overlap `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s established bounds on the independent Xenium breast-cancer dataset, after a subsampling bug was found and fixed.",
        "overall_evidence_grade": "strong",
    },
    {
        "claim_id": "q2_variance_partition_confirmatory",
        "discovery_evidence": "`13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`: nested variance-partition model (patient=29.3%, identity=20.2%, context=50.4%, bootstrap CIs).",
        "sensitivity_evidence": "`13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R`: clone-size confounder check, negligible shifts (<0.4pp).",
        "replicate_evidence": "`13_clone_ecology_confirmatory_models/04_run_leave_one_patient_out_stability.R`: leave-one-patient-out stability, 10/10 folds stable. `13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py`: segmentation-robustness check, engagement metrics concordant.",
        "external_validation_evidence": "Indirect only, via the underlying descriptors: `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` (module coherence transfers, GSE103322) and `16_external_validation_and_generalisation/03_validate_cell_state_signatures.R` (4/4 programmes significant, GSE139324) validate the T-cell state signatures this score is built from, not the variance-partition result itself.",
        "overall_evidence_grade": "strong",
    },
    {
        "claim_id": "q2_discrete_vs_continuous_structure_test",
        "discovery_evidence": "`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`: prespecified BIC + gap-statistic dual test -- continuous structure (both criteria required to agree, by design).",
        "sensitivity_evidence": "n/a -- the dual-criterion (BIC and gap) design is itself the sensitivity safeguard against a single statistic's own idiosyncrasies.",
        "replicate_evidence": "n/a -- a single cohort-wide structural test, not repeated per patient/replicate by design.",
        "external_validation_evidence": "`16_external_validation_and_generalisation/06_compare_with_source_paper_results.py`: diverges from the source paper's Figure 2 (14 discrete transcriptional clusters) -- an unreconciled scope/modality difference (targeted spatial panel vs. full-transcriptome scRNA).",
        "overall_evidence_grade": "moderate",
    },
    {
        "claim_id": "q3_barrier_topology_confirmatory",
        "discovery_evidence": "`14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`: nested mixed model, suppressive-myeloid barrier fraction significantly predicts lower engagement (estimate=-0.343, LRT p=0.0069, n=152).",
        "sensitivity_evidence": "n/a within Spatial Interactions and Barriers itself; the model already adjusts for state and niche composition as covariates.",
        "replicate_evidence": "`13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py` addendum: segmentation-robustness check on the same two barrier covariates, 2/2 same-sign in vs. outside the resegmented subset.",
        "external_validation_evidence": "`14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R`: direction-concordant with Grout et al. 2022's Fig 6B correlation (r=-0.48), but this cohort's raw correlation is notably weaker (r=-0.079, ratio 0.165) -- a partial concordance, not a full quantitative match.",
        "overall_evidence_grade": "moderate",
    },
    {
        "claim_id": "q3_literature_benchmark",
        "discovery_evidence": "n/a -- this claim IS itself external-validation evidence for q3_barrier_topology_confirmatory, not a discovery claim of its own.",
        "sensitivity_evidence": "n/a.",
        "replicate_evidence": "n/a.",
        "external_validation_evidence": "`14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R`: citations and reported statistics from Grout et al. 2022 (Cancer Discovery) and Hwang et al. 2022 (Nature Genetics).",
        "overall_evidence_grade": "supporting",
    },
    {
        "claim_id": "hpv_primary_contrast",
        "discovery_evidence": "`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`: validated n=4 vs n=4 contrast (2 patients excluded as clinically/molecularly discordant, `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`). `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`: 0/25 tests BH-significant.",
        "sensitivity_evidence": "`15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R`: LOPO (0/8 direction flips for the 2 stress-tested leads), bootstrap CIs, exhaustive permutation test -- all consistent, all short of formal significance.",
        "replicate_evidence": "n/a -- most contrast patients contribute a single primary section; the 2-section patients are already equal-weighted within `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`'s own patient-level aggregation.",
        "external_validation_evidence": "n/a -- HPV status is this claim's own exposure variable, not independently externally validated.",
        "overall_evidence_grade": "weak_exploratory",
    },
    {
        "claim_id": "null_model_calibration_suite",
        "discovery_evidence": "`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`: Type I error / power calibration on `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s synthetic ground truth, 3 null models x 7 effect sizes x 10 synthetic patients.",
        "sensitivity_evidence": "n/a -- is the sensitivity/calibration check itself.",
        "replicate_evidence": "10 synthetic-patient replicates.",
        "external_validation_evidence": "`16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py`: clean replication on an independent dataset after a subsampling bug fix (see q1_framework_generalisation, the same underlying evidence).",
        "overall_evidence_grade": "strong",
    },
    {
        "claim_id": "tcr_false_positive_rate",
        "discovery_evidence": "`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`: empirical false-positive-rate estimate, median 0.216, from three controls (off-patient probes, non-T cells, spatial autocorrelation of detection status).",
        "sensitivity_evidence": "The off-patient and non-T-cell rate-based controls are correlated across probes (Pearson r=0.528, Spearman rho=0.641, n=105 ), not fully independent; the spatial-autocorrelation control is a structurally different measurement (spatial clustering, not a detection rate) and is not implicated by this correlation. Three controls of two different kinds converging on the same estimate range remains informative, but is a weaker sensitivity argument than three mutually independent methods would be.",
        "replicate_evidence": "n/a.",
        "external_validation_evidence": "n/a.",
        "overall_evidence_grade": "strong",
    },
    {
        "claim_id": "replicate_concordance",
        "discovery_evidence": "`04_quality_control/08_assess_replicate_concordance.R`: 7 technical replicate pairs, 0/7 formally flagged discordant; P13 a convergent soft outlier across 3 independent metrics.",
        "sensitivity_evidence": "n/a -- is itself a replicate/sensitivity check.",
        "replicate_evidence": "The 7 replicate pairs themselves.",
        "external_validation_evidence": "`16_external_validation_and_generalisation/06_compare_with_source_paper_results.py`: diverges from the source paper's same-metric claim (they report 1/7 discordant under CDR3-probe Jaccard; this cohort's re-derivation finds 0/7) -- an investigated and explained divergence (binary Jaccard cannot see the quantitative-representation effect the source paper's text attributes the discordance to).",
        "overall_evidence_grade": "moderate",
    },
    {
        "claim_id": "segmentation_robustness_check",
        "discovery_evidence": "n/a -- this claim is itself a robustness/sensitivity check for q2_variance_partition_confirmatory and q3_barrier_topology_confirmatory, not a discovery claim of its own.",
        "sensitivity_evidence": "`13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py` + addendum: engagement and barrier metrics compared in vs. outside `04_quality_control/05_resegment_reference_subset.py`'s 3 resegmented sections, all metrics same-sign.",
        "replicate_evidence": "`04_quality_control/05_resegment_reference_subset.py`'s 3 resegmented sections (an independent transcript-reassignment pipeline).",
        "external_validation_evidence": "n/a.",
        "overall_evidence_grade": "supporting",
    },
    {
        "claim_id": "external_checkpoint_directional_consistency",
        "discovery_evidence": "`12_external_checkpoint_validation/02_quantify_directional_consistency.py`: pairwise sign agreement and rank-shift analysis between this cohort's T-cell state proportions and GSE103322's -- 80% pairwise directional agreement excluding the 2 discrepant states (Cycling, Ambiguous).",
        "sensitivity_evidence": "`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`'s module-coherence test (a different statistical question -- is the underlying gene module coherent, not just its abundance) corroborates that the Cycling discrepancy is an abundance effect, not a broken signature.",
        "replicate_evidence": "n/a.",
        "external_validation_evidence": "GSE103322 (External Checkpoint Validation, the claim's own external reference) + `16_external_validation_and_generalisation/03_validate_cell_state_signatures.R` (GSE139324, a second, independent reference, 4/4 programmes significant) -- two independent external checks, not one.",
        "overall_evidence_grade": "moderate",
    },
]


def validate_evidence_rows(rows: list[dict]) -> None:
    """Structural check: every row has all required non-empty fields and
    a valid `overall_evidence_grade`."""
    required_fields = [
        "claim_id",
        "discovery_evidence",
        "sensitivity_evidence",
        "replicate_evidence",
        "external_validation_evidence",
        "overall_evidence_grade",
    ]
    for row in rows:
        for field in required_fields:
            if not row.get(field):
                raise PipelineError(
                    f"Evidence row '{row.get('claim_id', '<unknown>')}' is missing a non-empty '{field}'."
                )
        if row["overall_evidence_grade"] not in VALID_GRADES:
            raise PipelineError(
                f"Evidence row '{row['claim_id']}' has an invalid overall_evidence_grade '{row['overall_evidence_grade']}'."
            )


def build_claim_evidence_matrix(project_root: Path) -> dict:
    registry_path = project_root / "governance" / "analysis_registry.tsv"
    output_path = project_root / "results" / "claim_evidence_matrix.tsv"

    if not registry_path.exists():
        raise PipelineError(f"'{registry_path}' not found.")

    registry = pd.read_csv(registry_path, sep="\t")
    registered_ids = set(registry["analysis_id"])
    covered_ids = {row["claim_id"] for row in CLAIM_EVIDENCE_MATRIX}
    missing = registered_ids - covered_ids
    if missing:
        raise PipelineError(
            f"Registered analysis/analyses {missing} in '{registry_path}' have no evidence row."
        )

    validate_evidence_rows(CLAIM_EVIDENCE_MATRIX)

    result = pd.DataFrame(CLAIM_EVIDENCE_MATRIX)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    grade_counts = result["overall_evidence_grade"].value_counts().to_dict()

    return {
        "n_claims": len(result),
        "grade_counts": grade_counts,
        "output_path": str(output_path),
    }
