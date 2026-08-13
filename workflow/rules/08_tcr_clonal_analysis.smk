"""TCR Clonal Analysis (08_tcr_clonal_analysis): TCR probe decoding, clone assignment, validation and ascertainment audit."""

rule phase08_00_identify_tcr_cdr3_probe_features:
    """08.00 -- Separates patient-specific CDR3 probes from conventional T-cell genes using feature metadata and naming rules. Primary output: metadata/tcr_probe_registry.tsv."""
    input: "results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/00_identify_tcr_cdr3_probe_features.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py"


rule phase08_01_map_tcr_probes_to_patients:
    """08.01 -- Ensures patient-specific probes are only evaluated in their intended specimens and detects leakage or naming conflicts. Primary output: reports/tcr/patient_probe_audit.tsv."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/00_identify_tcr_cdr3_probe_features.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/01_map_tcr_probes_to_patients.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py"


rule phase08_02_document_clone_ascertainment:
    """08.02 -- Records how each probed clonotype was selected relative to each patient's full repertoire and publishes this as an explicit boundary condition on every downstream generalisability claim. Primary output: metadata/clone_ascertainment.tsv."""
    input: "results/logs/08_tcr_clonal_analysis/.sentinels/01_map_tcr_probes_to_patients.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/02_document_clone_ascertainment.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/02_document_clone_ascertainment.py"


rule phase08_03_call_cell_level_tcr_detections:
    """08.03 -- Calls clone detections from transcript counts with explicit thresholds, controls and multi-probe ambiguity handling. Primary output: data/derived/tcr_cell_calls.parquet."""
    input: "results/logs/08_tcr_clonal_analysis/.sentinels/02_document_clone_ascertainment.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/03_call_cell_level_tcr_detections.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py"


rule phase08_04_estimate_false_positive_tcr_calls:
    """08.04 -- Uses off-patient probes, non-T cells and spatial permutation as empirical negative controls. Primary output: reports/tcr/false_positive_model.pdf."""
    input: "results/logs/08_tcr_clonal_analysis/.sentinels/02_document_clone_ascertainment.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/04_estimate_false_positive_tcr_calls.done")
    shell: "Rscript scripts/08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R"


rule phase08_05_screen_cdr3_cross_patient_similarity:
    """08.05 -- Screens probed CDR3 sequences for cross-patient similarity (public/quasi-public TCR motifs) that could cause probe cross-reactivity. Primary output: reports/tcr/cdr3_similarity_screen.tsv."""
    input: "results/logs/08_tcr_clonal_analysis/.sentinels/00_identify_tcr_cdr3_probe_features.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/05_screen_cdr3_cross_patient_similarity.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/05_screen_cdr3_cross_patient_similarity.py"


rule phase08_06_resolve_multiclonal_and_ambiguous_cells:
    """08.06 -- Classifies singlet, probable multiplet, low-confidence and unassigned TCR calls. Primary output: data/derived/tcr_resolved_calls.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/03_call_cell_level_tcr_detections.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/04_estimate_false_positive_tcr_calls.done",
        "results/logs/04_quality_control/.sentinels/04_estimate_transcript_spillover.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/06_resolve_multiclonal_and_ambiguous_cells.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py"


rule phase08_07_build_clone_metadata_table:
    """08.07 -- Summarises clone size, patient, section support, phenotype composition and replicate recurrence. Primary output: data/derived/clone_metadata.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/06_resolve_multiclonal_and_ambiguous_cells.done",
        "results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/07_build_clone_metadata_table.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/07_build_clone_metadata_table.py"


rule phase08_08_generate_tcr_release_report:
    """08.08 -- Freezes high-confidence clone definitions for primary analysis and documents excluded/ambiguous calls. Primary output: data/releases/v1_tcr_calls/."""
    input: "results/logs/08_tcr_clonal_analysis/.sentinels/07_build_clone_metadata_table.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/08_generate_tcr_release_report.py"


rule phase08_09_validate_probe_clones_against_paired_vdj_ground_truth:
    """08.09 -- Post-hoc follow-up, not part of the original 18-stage plan (see docs/analysis_amendments.md): validates Xenium CDR3-probe detections against independent, paired scTCR-seq VDJ ground truth from GSE287301 (same 28-patient cohort) -- directly superseding an earlier 'no VDJ data available' finding. Primary output: data/derived/probe_vdj_ground_truth_comparison.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/01_map_tcr_probes_to_patients.done",
        "results/logs/06_cell_type_annotation/.sentinels/08_acquire_companion_scrna_and_vdj_reference.done"
    output: touch("results/logs/08_tcr_clonal_analysis/.sentinels/09_validate_probe_clones_against_paired_vdj_ground_truth.done")
    shell: "python3 scripts/08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_vdj_ground_truth.py"
