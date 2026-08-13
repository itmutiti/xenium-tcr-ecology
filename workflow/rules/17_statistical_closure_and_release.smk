"""Statistical Closure and Release (17_statistical_closure_and_release): Statistical closure, figures, manuscript and reproducible release."""

rule phase17_00_freeze_primary_results:
    """17.00 -- Locks the prespecified primary results (Q1-Q3, plus the single HPV contrast set) before exploratory extensions. Primary output: data/releases/final_primary/."""
    input:
        "results/logs/16_external_validation_and_generalisation/.sentinels/05_validate_framework_on_independent_dataset.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/01_fit_variance_partition_models.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/05_test_discrete_vs_continuous_structure.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/03_model_barrier_topology_by_structure.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/03_compare_cellular_composition_patient_level.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/04_compare_ecosystem_and_clone_structure.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/05_run_small_sample_robustness_checks.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/00_freeze_primary_results.done")
    shell: "python3 scripts/17_statistical_closure_and_release/00_freeze_primary_results.py"


rule phase17_01_control_multiplicity_and_report_effects:
    """17.01 -- Applies a single, project-wide FDR/gatekeeping strategy across the prespecified primary family. Primary output: results/statistical_summary.tsv."""
    input: "results/logs/17_statistical_closure_and_release/.sentinels/00_freeze_primary_results.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/01_control_multiplicity_and_report_effects.done")
    shell: "Rscript scripts/17_statistical_closure_and_release/01_control_multiplicity_and_report_effects.R"


rule phase17_11_build_redesigned_manuscript_figures:
    """17.11 -- Builds 4 composite manuscript figures from already-validated data (framework generalisation on 3 tumour types, variance partition with its sensitivity check, barrier topology with its covariate-ablation check, consolidated HPV); see docs/analysis_amendments.md. Must run before 17.02, which assembles these into figures/main/. Primary output: reports/manuscript_figures/."""
    input:
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/08_run_calibration_suite_on_synthetic_data.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/05_validate_framework_on_independent_dataset.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/09_validate_framework_on_second_cancer_type.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/01_fit_variance_partition_models.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/07_test_structure_sensitivity_excluding_cycling.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/03_model_barrier_topology_by_structure.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/06_benchmark_against_published_barrier_studies.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/07_ablate_covariates_for_barrier_effect.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/03_compare_cellular_composition_patient_level.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/04_compare_ecosystem_and_clone_structure.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/11_build_redesigned_manuscript_figures.done")
    shell: "python3 scripts/17_statistical_closure_and_release/11_build_redesigned_manuscript_figures.py"


rule phase17_02_generate_all_main_figures:
    """17.02 -- Generates manuscript figures from frozen data with no manual editing of quantitative content. Primary output: figures/main/."""
    input: "results/logs/17_statistical_closure_and_release/.sentinels/11_build_redesigned_manuscript_figures.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/02_generate_all_main_figures.done")
    shell: "python3 scripts/17_statistical_closure_and_release/02_generate_all_main_figures.py"


rule phase17_03_generate_all_supplementary_figures:
    """17.03 -- Creates complete diagnostics, sensitivity analyses, segmentation-robustness comparisons and non-selected clone maps. Primary output: figures/supplementary/."""
    input:
        "results/logs/.checkpoints/full_validation_release.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/06_generate_clone_atlas.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/03_generate_all_supplementary_figures.done")
    shell: "python3 scripts/17_statistical_closure_and_release/03_generate_all_supplementary_figures.py"


rule phase17_04_generate_results_tables:
    """17.04 -- Creates publication-ready sample, QC, model and validation tables. Primary output: tables/."""
    input:
        "results/logs/17_statistical_closure_and_release/.sentinels/01_control_multiplicity_and_report_effects.done",
        "results/logs/04_quality_control/.sentinels/08_assess_replicate_concordance.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/00_define_validation_claims.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/06_compare_with_source_paper_results.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/07_generate_evidence_matrix.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/06_prepare_hpv_claim_strength_table.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/04_generate_results_tables.done")
    shell: "python3 scripts/17_statistical_closure_and_release/04_generate_results_tables.py"


rule phase17_05_render_reproducible_manuscript:
    """17.05 -- Builds the manuscript with embedded software and data provenance from `05_render_reproducible_manuscript.qmd`'s source content. Primary output: manuscript/manuscript.pdf. Rendered by `05b_render_manuscript_pdf.py` via a pure-Python markdown + xhtml2pdf path. See src/xenium_tcr_ecology/release/manuscript.py's module docstring."""
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/05_render_reproducible_manuscript.done")
    shell: "python3 scripts/17_statistical_closure_and_release/05b_render_manuscript_pdf.py"


rule phase17_06_build_public_data_release:
    """17.06 -- Exports non-sensitive derived data, schemas and lightweight examples with licences. Primary output: release/data/. Runs only under the publication_release target, not the default one; see workflow/rules/_checkpoints.smk's publication_release target."""
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/06_build_public_data_release.done")
    shell: "python3 scripts/17_statistical_closure_and_release/06_build_public_data_release.py"


rule phase17_07_build_documented_software_package:
    """17.07 -- Packages reusable functions, API documentation, tutorials and tests as an installable, versioned package. Primary output: release/software/."""
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/07_build_documented_software_package.done")
    shell: "python3 scripts/17_statistical_closure_and_release/07_build_documented_software_package.py"


rule phase17_08_run_end_to_end_reproduction_test:
    """17.08 -- Executes the full pipeline inside the container image, on a clean environment, first on a test subset and then the full dataset. Primary output: reports/reproduction_test/.
    NOTE: this smoke test's own script docstring states it runs "as part of
    computational_analysis_complete" -- i.e. late, after the derived outputs
    its representative script subset needs already exist (see
    docs/execution_manual/EXECUTION_MANUAL.md's "Smoke test vs. full
    reproduction" section). input: below reflects that: it depends on
    results tables being assembled, the last regular computational step
    before closure."""
    input: "results/logs/17_statistical_closure_and_release/.sentinels/04_generate_results_tables.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/08_run_end_to_end_reproduction_test.done")
    shell: "bash scripts/17_statistical_closure_and_release/08_run_end_to_end_reproduction_test.sh"


rule phase17_09_run_null_model_calibration_regression:
    """17.09 -- Re-runs the Spatial Graph Construction and Calibration simulation-calibration suite against the frozen release code, guarding against calibration drift. Primary output: reports/release/calibration_regression.tsv."""
    input: "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/08_run_calibration_suite_on_synthetic_data.done"
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/09_run_null_model_calibration_regression.done")
    shell: "python3 scripts/17_statistical_closure_and_release/09_run_null_model_calibration_regression.py"


rule phase17_10_archive_code_and_create_doi:
    """17.10 -- Creates a tagged release, citation file, software bill of materials, container image digest and archive-ready bundle. Primary output: release/archive_bundle/. Creates a permanent, citable DOI; runs only under the publication_release target, not the default one; see workflow/rules/_checkpoints.smk's publication_release target."""
    output: touch("results/logs/17_statistical_closure_and_release/.sentinels/10_archive_code_and_create_doi.done")
    shell: "bash scripts/17_statistical_closure_and_release/10_archive_code_and_create_doi.sh"
