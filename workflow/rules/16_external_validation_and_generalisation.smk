"""External Validation and Generalisation (16_external_validation_and_generalisation): Full external validation, triangulation and generalisability."""

rule phase16_00_define_validation_claims:
    """16.00 -- Maps each major claim to a specific independent validation dataset and success criterion, predeclared before analysis. Primary output: governance/validation_plan.tsv."""
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/00_define_validation_claims.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/00_define_validation_claims.py"


rule phase16_01_acquire_independent_spatial_dataset:
    """16.01 -- Seeks at least one independent spatial dataset with clonal or lineage information to test whether the framework generalises. Primary output: data/external/spatial/."""
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/01_acquire_independent_spatial_dataset.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py"


rule phase16_02_acquire_hnscc_scrna_references:
    """16.02 -- Downloads selected public HNSCC scRNA-seq datasets with licences, metadata and checksums. Primary output: data/external/scrna/."""
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/02_acquire_hnscc_scrna_references.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py"


rule phase16_03_validate_cell_state_signatures:
    """16.03 -- Tests whether targeted-panel programs correspond to full-transcriptome T-cell and stromal states. Primary output: reports/validation/signature_validation.pdf."""
    input: "results/logs/16_external_validation_and_generalisation/.sentinels/02_acquire_hnscc_scrna_references.done"
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/03_validate_cell_state_signatures.done")
    shell: "Rscript scripts/16_external_validation_and_generalisation/03_validate_cell_state_signatures.R"


rule phase16_04_validate_ecosystem_signatures_in_bulk:
    """16.04 -- Projects ecosystem-derived signatures to TCGA/CPTAC-like cohorts and cautiously tests clinical/immune associations. Primary output: reports/validation/bulk_projection.pdf."""
    input: "results/logs/16_external_validation_and_generalisation/.sentinels/02_acquire_hnscc_scrna_references.done"
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/04_validate_ecosystem_signatures_in_bulk.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/04_validate_ecosystem_signatures_in_bulk.py"


rule phase16_05_validate_framework_on_independent_dataset:
    """16.05 -- Applies the full software framework end-to-end to the independent spatial dataset - the strongest single test of Q1. Primary output: reports/validation/framework_generalisation.pdf."""
    input:
        "results/logs/16_external_validation_and_generalisation/.sentinels/01_acquire_independent_spatial_dataset.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/02_calibrate_graph_parameters.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/08_run_calibration_suite_on_synthetic_data.done"
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/05_validate_framework_on_independent_dataset.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py"


rule phase16_06_compare_with_source_paper_results:
    """16.06 -- Reproduces selected published McCord et al. results as a benchmark, then identifies which findings in this thesis are new. Primary output: reports/validation/source_paper_comparison.tsv."""
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/06_compare_with_source_paper_results.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/06_compare_with_source_paper_results.py"


rule phase16_07_generate_evidence_matrix:
    """16.07 -- Links every manuscript claim to discovery, sensitivity, replicate and external-validation evidence. Primary output: results/claim_evidence_matrix.tsv."""
    input: "results/logs/01_project_setup_and_governance/.sentinels/06_create_analysis_registry.done"
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/07_generate_evidence_matrix.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/07_generate_evidence_matrix.py"


rule phase16_08_acquire_second_independent_spatial_dataset:
    """16.08 -- Post-hoc follow-up, not part of the original 18-stage plan (see docs/analysis_amendments.md): verifies the second, independent Xenium colorectal-cancer dataset (de Oliveira et al. 2025) used to strengthen q1_framework_generalisation. Primary output: data/external/spatial/Xenium_Oliveira_ColorectalCancer_P1/."""
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/08_acquire_second_independent_spatial_dataset.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/08_acquire_second_independent_spatial_dataset.py"


rule phase16_09_validate_framework_on_second_cancer_type:
    """16.09 -- Post-hoc follow-up (see docs/analysis_amendments.md): applies the calibrated null-model framework to the second, independent colorectal-cancer Xenium dataset, testing whether calibration holds across more than one tissue type.Primary output: reports/validation/framework_generalisation_second_dataset.pdf."""
    input:
        "results/logs/16_external_validation_and_generalisation/.sentinels/08_acquire_second_independent_spatial_dataset.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/02_calibrate_graph_parameters.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/08_run_calibration_suite_on_synthetic_data.done"
    output: touch("results/logs/16_external_validation_and_generalisation/.sentinels/09_validate_framework_on_second_cancer_type.done")
    shell: "python3 scripts/16_external_validation_and_generalisation/09_validate_framework_on_second_cancer_type.py"
