"""Preprocessing and Normalisation (05_preprocessing_and_normalisation): Expression normalisation, feature engineering and technical variation."""

rule phase05_00_separate_gene_and_control_features:
    """05.00 -- Separates biological genes, negative controls, HPV probes and TCR/CDR3 probes into explicit feature classes. Primary output: metadata/feature_annotation.tsv."""
    input: "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done")
    shell: "python3 scripts/05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py"


rule phase05_01_construct_analysis_count_layers:
    """05.01 -- Preserves raw counts and creates normalised, variance-stabilised and binary-detection layers where appropriate. Primary output: layers in h5ad."""
    input:
        "results/logs/04_quality_control/.sentinels/07_apply_qc_filters_with_audit_trail.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done"
    output: touch("results/logs/05_preprocessing_and_normalisation/.sentinels/01_construct_analysis_count_layers.done")
    shell: "python3 scripts/05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py"


rule phase05_02_evaluate_normalisation_strategies:
    """05.02 -- Compares log-normalisation, Pearson residuals and detection-based approaches using replicate stability and negative-control probe behaviour. Primary output: reports/preprocess/normalisation_benchmark.pdf."""
    input: "results/logs/05_preprocessing_and_normalisation/.sentinels/01_construct_analysis_count_layers.done"
    output: touch("results/logs/05_preprocessing_and_normalisation/.sentinels/02_evaluate_normalisation_strategies.done")
    shell: "Rscript scripts/05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R"


rule phase05_03_calculate_program_scores:
    """05.03 -- Computes curated cytotoxicity, exhaustion, activation, interferon, proliferation, stress, EMT and antigen-presentation scores. Primary output: data/derived/program_scores.parquet.
    Needs 02_evaluate_normalisation_strategies.R's uns['primary_normalization_layer']
    choice on the shared analysis_ready.h5ad -- found missing during the second Vast.ai clean-room run: both
    05.02 and 05.03 depended only on 05.01, so parallel scheduling could (and
    did) run 05.03 before 05.02 had written that attribute."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/01_construct_analysis_count_layers.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/02_evaluate_normalisation_strategies.done"
    output: touch("results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done")
    shell: "python3 scripts/05_preprocessing_and_normalisation/03_calculate_program_scores.py"


rule phase05_04_model_technical_covariates:
    """05.04 -- Estimates contributions of run, section, depth and control burden without erasing patient biology; run is modelled explicitly as nested within patient. Primary output: reports/preprocess/variance_partition.pdf."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done",
        "results/logs/04_quality_control/.sentinels/00_compute_cell_level_qc_metrics.done"
    output: touch("results/logs/05_preprocessing_and_normalisation/.sentinels/04_model_technical_covariates.done")
    shell: "Rscript scripts/05_preprocessing_and_normalisation/04_model_technical_covariates.R"


rule phase05_05_create_primary_analysis_matrix:
    """05.05 -- Freezes the matrix and covariates used by downstream analyses and records its exact hash. Primary output: data/releases/v1_primary_analysis/."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/01_construct_analysis_count_layers.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done"
    output: touch("results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done")
    shell: "python3 scripts/05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py"
