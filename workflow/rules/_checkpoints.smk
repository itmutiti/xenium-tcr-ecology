"""Named checkpoints for the reviewer-facing computational DAG (raw data ->
final biological results)."""

rule project_ready:
    """Deliberately does NOT require 04_build_container_images.done.
    `scripts/01_project_setup_and_governance/04_build_container_images.sh`
    hard-fails if a Docker daemon isn't reachable, with no fallback --
    confirmed absent on both Vast.ai clean-room instances used to date
    (unprivileged containers, no CAP_SYS_ADMIN). Nothing downstream reads
    the built image or `containers/release.sif`: checked by code read --
    `src/xenium_tcr_ecology/release/reproduction_test.py`'s own Docker
    mention is historical narrative describing a past manual run, not a
    live dependency (its actual smoke test invokes `sys.executable`
    directly, never `docker run`). Requiring it here would block the
    entire computational DAG on packaging infrastructure that nothing
    scientific depends on -- the same class of spurious blocking
    dependency already fixed for phase12_04 and phase01_00's original
    flag-passthrough gap. `phase01_04_build_container_images` remains a
    fully real, unmodified, independently-runnable rule (not deleted):
    run it by name on a host with Docker support when a
    container image is actually needed."""
    input:
        "results/logs/01_project_setup_and_governance/.sentinels/00_validate_project_scope.done",
        "results/logs/01_project_setup_and_governance/.sentinels/01_build_sample_manifest.done",
        "results/logs/01_project_setup_and_governance/.sentinels/02_generate_data_dictionary.done",
        "results/logs/01_project_setup_and_governance/.sentinels/03_lock_software_environments.done",
        "results/logs/01_project_setup_and_governance/.sentinels/05_initialise_reproducible_workflow.done",
        "results/logs/01_project_setup_and_governance/.sentinels/06_create_analysis_registry.done"
    output: touch("results/logs/.checkpoints/project_ready.done")


rule raw_data_frozen:
    input:
        "results/logs/.checkpoints/project_ready.done",
        "results/logs/02_raw_data_ingestion/.sentinels/00_query_geo_accession.done",
        "results/logs/02_raw_data_ingestion/.sentinels/01_download_geo_raw_archive.done",
        "results/logs/02_raw_data_ingestion/.sentinels/02_verify_archive_checksums.done",
        "results/logs/02_raw_data_ingestion/.sentinels/03_extract_archive_safely.done",
        "results/logs/02_raw_data_ingestion/.sentinels/04_inventory_xenium_files.done",
        "results/logs/02_raw_data_ingestion/.sentinels/05_standardise_sample_directory_layout.done",
        "results/logs/02_raw_data_ingestion/.sentinels/06_freeze_raw_data_snapshot.done"
    output: touch("results/logs/.checkpoints/raw_data_frozen.done")


rule canonical_objects:
    input:
        "results/logs/.checkpoints/raw_data_frozen.done",
        "results/logs/03_spatialdata_import/.sentinels/00_detect_xenium_format_version.done",
        "results/logs/03_spatialdata_import/.sentinels/01_import_each_section_to_spatialdata.done",
        "results/logs/03_spatialdata_import/.sentinels/02_validate_coordinate_systems.done",
        "results/logs/03_spatialdata_import/.sentinels/03_create_anndata_expression_objects.done",
        "results/logs/03_spatialdata_import/.sentinels/04_attach_clinical_and_technical_metadata.done",
        "results/logs/03_spatialdata_import/.sentinels/05_build_combined_analysis_object.done",
        "results/logs/03_spatialdata_import/.sentinels/06_export_r_interoperability_objects.done"
    output: touch("results/logs/.checkpoints/canonical_objects.done")


rule qc_and_integrity_release:
    """Gates entry to Preprocessing and Normalisation / Cell Type Annotation (annotation_release) on
    Quality Control's core threshold/filter/release pipeline (4.00-4.03, 4.06-4.09)
    only -- NOT on 4.04 (transcript spillover) or 4.05 (resegmentation
    robustness), which are deliberately deferred until after Cell Type Annotation
    produces real cell-type labels: neither gates the QC threshold/filter/release
    pipeline, since 04.06-04.09 depend only on the already-computed
    04.00-04.03 metrics, not on spillover-risk or the resegmented subset.
    Requiring 4.04/4.05 here would be a cyclic dependency (they need Cell Type Annotation,
    which needs this checkpoint) and was never a requirement for what
    this checkpoint actually gates. 4.04/4.05's real dependents
    (phase08_06_resolve_multiclonal_and_ambiguous_cells,
    phase11_02_quantify_clone_tumour_engagement,
    phase13_05_test_segmentation_robustness) depend on their
    sentinels directly instead, once those phases are reached."""
    input:
        "results/logs/.checkpoints/canonical_objects.done",
        "results/logs/04_quality_control/.sentinels/00_compute_cell_level_qc_metrics.done",
        "results/logs/04_quality_control/.sentinels/01_compute_transcript_level_qc_metrics.done",
        "results/logs/04_quality_control/.sentinels/02_detect_spatial_qc_artifacts.done",
        "results/logs/04_quality_control/.sentinels/03_assess_segmentation_quality.done",
        "results/logs/04_quality_control/.sentinels/06_define_qc_thresholds_hierarchically.done",
        "results/logs/04_quality_control/.sentinels/07_apply_qc_filters_with_audit_trail.done",
        "results/logs/04_quality_control/.sentinels/08_assess_replicate_concordance.done",
        "results/logs/04_quality_control/.sentinels/09_generate_qc_release_report.done"
    output: touch("results/logs/.checkpoints/qc_and_integrity_release.done")


rule annotation_release:
    input:
        "results/logs/.checkpoints/qc_and_integrity_release.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/01_construct_analysis_count_layers.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/02_evaluate_normalisation_strategies.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/04_model_technical_covariates.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/00_compile_marker_and_reference_registry.done",
        "results/logs/06_cell_type_annotation/.sentinels/01_cluster_within_patient_and_jointly.done",
        "results/logs/06_cell_type_annotation/.sentinels/02_score_major_lineages.done",
        "results/logs/06_cell_type_annotation/.sentinels/03_map_external_scrna_reference.done",
        "results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done",
        "results/logs/06_cell_type_annotation/.sentinels/05_resolve_myeloid_and_stromal_substates.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/06_cell_type_annotation/.sentinels/07_blinded_annotation_review.done"
    output: touch("results/logs/.checkpoints/annotation_release.done")


rule tumour_and_tcr_release:
    input:
        "results/logs/.checkpoints/annotation_release.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/00_subset_and_recluster_epithelial_cells.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/01_score_malignancy_and_normal_epithelium.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/02_cross_validate_against_morphology.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/03_infer_cnv_appendix_only.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/04_construct_tumour_region_masks.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/05_extract_tumour_boundaries.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/06_validate_boundaries_against_morphology.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/07_define_invasive_front_and_compartments.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/00_identify_tcr_cdr3_probe_features.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/01_map_tcr_probes_to_patients.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/02_document_clone_ascertainment.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/03_call_cell_level_tcr_detections.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/04_estimate_false_positive_tcr_calls.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/05_screen_cdr3_cross_patient_similarity.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/06_resolve_multiclonal_and_ambiguous_cells.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/07_build_clone_metadata_table.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/09_validate_probe_clones_against_paired_vdj_ground_truth.done"
    output: touch("results/logs/.checkpoints/tumour_and_tcr_release.done")


rule graph_and_calibration_release:
    input:
        "results/logs/.checkpoints/tumour_and_tcr_release.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/00_generate_candidate_spatial_graphs.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/01_prune_graphs_for_tissue_gaps.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/02_calibrate_graph_parameters.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/04_construct_tumour_tcell_bipartite_graph.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/05_construct_clone_induced_subgraphs.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/06_run_graph_sensitivity_grid.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/07_generate_synthetic_ground_truth_patterns.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/08_run_calibration_suite_on_synthetic_data.done"
    output: touch("results/logs/.checkpoints/graph_and_calibration_release.done")


rule niche_release:
    input:
        "results/logs/.checkpoints/graph_and_calibration_release.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/00_compute_cell_type_neighbourhood_enrichment.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/02_discover_neighbourhood_archetypes.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/03_segment_tissue_domains.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/04_annotate_ecosystems_with_blinded_rules.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/05_quantify_ecosystem_abundance_and_topology.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/06_test_patient_recurrence.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/07_leave_one_patient_out_niche_stability.done"
    output: touch("results/logs/.checkpoints/niche_release.done")


rule clone_descriptor_release:
    input:
        "results/logs/.checkpoints/niche_release.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/00_compute_clone_spatial_descriptors_rarefied.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/01_compute_clone_cell_state_composition.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/03_quantify_clone_apc_support.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/05_test_discrete_vs_continuous_structure.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/06_discover_provisional_structure.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/07_freeze_provisional_taxonomy_version.done"
    output: touch("results/logs/.checkpoints/clone_descriptor_release.done")


rule external_checkpoint:
    """12.03's predeclared, fully automated decision rule (freeze iff all
    programs transfer AND overall Spearman rho > 0; revise otherwise) is
    the taxonomy-freeze decision -- there is no separate attribution or
    approval step in this DAG. Its output, governance/freeze_decision.tsv
    (decision=freeze), is read directly by
    src/xenium_tcr_ecology/clone_ecology/taxonomy_release.py; this
    checkpoint requires nothing beyond 12.03 completing."""
    input:
        "results/logs/.checkpoints/clone_descriptor_release.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/00_project_provisional_signatures_to_bulk_reference.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/01_test_transcriptional_program_transfer.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/02_quantify_directional_consistency.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/03_decide_freeze_or_revise.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/05_rescore_cycling_state_with_primary_method.done"
    output: touch("results/logs/.checkpoints/external_checkpoint.done")


rule clone_model_release:
    input:
        "results/logs/.checkpoints/external_checkpoint.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/01_fit_variance_partition_models.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/02_fit_hierarchical_clone_models.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/03_test_clone_size_as_confounder.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/04_run_leave_one_patient_out_stability.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/05_test_segmentation_robustness.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/06_generate_clone_atlas.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/07_test_structure_sensitivity_excluding_cycling.done"
    output: touch("results/logs/.checkpoints/clone_model_release.done")


rule interaction_release:
    input:
        "results/logs/.checkpoints/clone_model_release.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/00_define_sender_receiver_pairs.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/01_filter_ligand_receptor_database_to_panel.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/02_compute_spatially_constrained_scores.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/03_model_barrier_topology_by_structure.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/04_analyse_barrier_pathways.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/05_prioritise_testable_interactions.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/06_benchmark_against_published_barrier_studies.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/07_ablate_covariates_for_barrier_effect.done"
    output: touch("results/logs/.checkpoints/interaction_release.done")


rule hpv_release:
    input:
        "results/logs/.checkpoints/interaction_release.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/00_validate_hpv_metadata_and_probe_signal.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/02_run_prospective_power_simulation.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/03_compare_cellular_composition_patient_level.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/04_compare_ecosystem_and_clone_structure.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/05_run_small_sample_robustness_checks.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/06_prepare_hpv_claim_strength_table.done"
    output: touch("results/logs/.checkpoints/hpv_release.done")


rule full_validation_release:
    input:
        "results/logs/.checkpoints/hpv_release.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/00_define_validation_claims.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/01_acquire_independent_spatial_dataset.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/02_acquire_hnscc_scrna_references.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/03_validate_cell_state_signatures.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/04_validate_ecosystem_signatures_in_bulk.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/05_validate_framework_on_independent_dataset.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/06_compare_with_source_paper_results.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/07_generate_evidence_matrix.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/08_acquire_second_independent_spatial_dataset.done",
        "results/logs/16_external_validation_and_generalisation/.sentinels/09_validate_framework_on_second_cancer_type.done"
    output: touch("results/logs/.checkpoints/full_validation_release.done")


rule computational_analysis_complete:
    """Final checkpoint of the reviewer-facing computational DAG (raw data
    -> final biological results, all statistical closure and manuscript-
    ready outputs). This is `rule all`'s default target -- running
    `snakemake --cores N` with no target stops here. Deliberately does
    not include 17_statistical_closure_and_release's publish-only steps (public data export,
    code archival/DOI creation, both external-facing actions) -- see
    `rule publication_release` below."""
    input:
        "results/logs/.checkpoints/full_validation_release.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/00_freeze_primary_results.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/01_control_multiplicity_and_report_effects.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/11_build_redesigned_manuscript_figures.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/02_generate_all_main_figures.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/03_generate_all_supplementary_figures.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/04_generate_results_tables.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/05_render_reproducible_manuscript.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/07_build_documented_software_package.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/08_run_end_to_end_reproduction_test.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/09_run_null_model_calibration_regression.done"
    output: touch("results/logs/.checkpoints/computational_analysis_complete.done")


rule publication_release:
    """Never run by default (`rule all` stops at
    `computational_analysis_complete` above) and never invoked by CI or
    any other automated process. Invoked by explicit name
    (`snakemake --cores N publication_release`). 17.06 exports data
    publicly and 17.10 creates a permanent, citable DOI archive
    (irreversible once registered)."""
    input:
        "results/logs/.checkpoints/computational_analysis_complete.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/06_build_public_data_release.done",
        "results/logs/17_statistical_closure_and_release/.sentinels/10_archive_code_and_create_doi.done"
    output: touch("results/logs/.checkpoints/publication_release.done")
