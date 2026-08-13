# Phase reference

Generated from `manifests/phase_registry.yaml` -- do not hand-edit.

## `01_project_setup_and_governance`

Research governance, computational environment and reproducibility

Deterministic: `True`

Scripts:

- `00_validate_project_scope.py`
- `01_build_sample_manifest.py`
- `02_generate_data_dictionary.py`
- `03_lock_software_environments.sh`
- `04_build_container_images.sh`
- `05_initialise_reproducible_workflow.py`
- `06_create_analysis_registry.py`

## `02_raw_data_ingestion`

Data acquisition, integrity verification and immutable staging

Deterministic: `True`

Scripts:

- `00_query_geo_accession.py`
- `01_download_geo_raw_archive.sh`
- `02_verify_archive_checksums.py`
- `03_extract_archive_safely.py`
- `04_inventory_xenium_files.py`
- `05_standardise_sample_directory_layout.py`
- `06_freeze_raw_data_snapshot.sh`

## `03_spatialdata_import`

Xenium import, coordinate harmonisation and canonical data objects

Deterministic: `True`

Scripts:

- `00_detect_xenium_format_version.py`
- `01_import_each_section_to_spatialdata.py`
- `02_validate_coordinate_systems.py`
- `03_create_anndata_expression_objects.py`
- `04_attach_clinical_and_technical_metadata.py`
- `05_build_combined_analysis_object.py`
- `06_export_r_interoperability_objects.py`

## `04_quality_control`

Multilevel quality control, principled exclusion and transcript-integrity robustness

Deterministic: `True`

Scripts:

- `00_compute_cell_level_qc_metrics.py`
- `01_compute_transcript_level_qc_metrics.py`
- `02_detect_spatial_qc_artifacts.py`
- `03_assess_segmentation_quality.py`
- `04_estimate_transcript_spillover.py`
- `05_resegment_reference_subset.py`
- `06_define_qc_thresholds_hierarchically.R`
- `07_apply_qc_filters_with_audit_trail.py`
- `08_assess_replicate_concordance.R`
- `09_generate_qc_release_report.py`

## `05_preprocessing_and_normalisation`

Expression normalisation, feature engineering and technical variation

Deterministic: `True`

Scripts:

- `00_separate_gene_and_control_features.py`
- `01_construct_analysis_count_layers.py`
- `02_evaluate_normalisation_strategies.R`
- `03_calculate_program_scores.py`
- `04_model_technical_covariates.R`
- `05_create_primary_analysis_matrix.py`
- `_02_compute_normalization_benchmark_metrics.py`

## `06_cell_type_annotation`

Hierarchical cell annotation with uncertainty

Deterministic: `True`

Scripts:

- `00_compile_marker_and_reference_registry.py`
- `01_cluster_within_patient_and_jointly.py`
- `02_score_major_lineages.py`
- `03_map_external_scrna_reference.py`
- `04_resolve_t_cell_substates.R`
- `05_resolve_myeloid_and_stromal_substates.R`
- `06_integrate_annotation_evidence.py`
- `07_blinded_annotation_review.py`
- `08_acquire_companion_scrna_and_vdj_reference.py`
- `_04_prepare_t_cell_substate_inputs.py`
- `_05_prepare_tme_substate_inputs.py`

## `07_tumour_epithelium_characterisation`

Malignant epithelial-state inference and tumour-boundary reconstruction

Deterministic: `True`

Scripts:

- `00_subset_and_recluster_epithelial_cells.py`
- `01_score_malignancy_and_normal_epithelium.py`
- `02_cross_validate_against_morphology.py`
- `03_infer_cnv_appendix_only.R`
- `04_construct_tumour_region_masks.py`
- `05_extract_tumour_boundaries.py`
- `06_validate_boundaries_against_morphology.py`
- `07_define_invasive_front_and_compartments.py`
- `_03_prepare_cnv_inputs.py`

## `08_tcr_clonal_analysis`

TCR probe decoding, clone assignment, validation and ascertainment audit

Deterministic: `True`

Scripts:

- `00_identify_tcr_cdr3_probe_features.py`
- `01_map_tcr_probes_to_patients.py`
- `02_document_clone_ascertainment.py`
- `03_call_cell_level_tcr_detections.py`
- `04_estimate_false_positive_tcr_calls.R`
- `05_screen_cdr3_cross_patient_similarity.py`
- `06_resolve_multiclonal_and_ambiguous_cells.py`
- `07_build_clone_metadata_table.py`
- `08_generate_tcr_release_report.py`
- `09_validate_probe_clones_against_paired_vdj_ground_truth.py`
- `_04_prepare_false_positive_inputs.py`

## `09_spatial_graph_construction_and_calibration`

Spatial graph construction, sensitivity framework and null-model calibration

Deterministic: `True`

Scripts:

- `00_generate_candidate_spatial_graphs.py`
- `01_prune_graphs_for_tissue_gaps.py`
- `02_calibrate_graph_parameters.py`
- `03_construct_primary_cell_graph.py`
- `04_construct_tumour_tcell_bipartite_graph.py`
- `05_construct_clone_induced_subgraphs.py`
- `06_run_graph_sensitivity_grid.py`
- `07_generate_synthetic_ground_truth_patterns.py`
- `08_run_calibration_suite_on_synthetic_data.py`

## `10_niche_and_ecosystem_discovery`

Discovery of spatial ecosystems and tissue domains

Deterministic: `True`

Scripts:

- `00_compute_cell_type_neighbourhood_enrichment.py`
- `01_compute_local_neighbourhood_compositions.py`
- `02_discover_neighbourhood_archetypes.R`
- `03_segment_tissue_domains.py`
- `04_annotate_ecosystems_with_blinded_rules.py`
- `05_quantify_ecosystem_abundance_and_topology.py`
- `06_test_patient_recurrence.R`
- `07_leave_one_patient_out_niche_stability.py`

## `11_clone_spatial_descriptors`

Clone-level ecological descriptors and structure test (provisional)

Deterministic: `True`

Scripts:

- `00_compute_clone_spatial_descriptors_rarefied.py`
- `01_compute_clone_cell_state_composition.py`
- `02_quantify_clone_tumour_engagement.py`
- `03_quantify_clone_apc_support.py`
- `04_quantify_stromal_and_myeloid_barriers.py`
- `05_test_discrete_vs_continuous_structure.R`
- `06_discover_provisional_structure.R`
- `07_freeze_provisional_taxonomy_version.py`

## `12_external_checkpoint_validation`

Early external sanity-check and taxonomy-freeze gate

Deterministic: `True`

Scripts:

- `00_project_provisional_signatures_to_bulk_reference.py`
- `01_test_transcriptional_program_transfer.py`
- `02_quantify_directional_consistency.py`
- `03_decide_freeze_or_revise.py`
- `05_rescore_cycling_state_with_primary_method.py`

## `13_clone_ecology_confirmatory_models`

Confirmatory variance-partitioning and clone ecological modelling

Deterministic: `True`

Scripts:

- `00_load_frozen_taxonomy_version.py`
- `01_fit_variance_partition_models.R`
- `02_fit_hierarchical_clone_models.R`
- `03_test_clone_size_as_confounder.R`
- `04_run_leave_one_patient_out_stability.R`
- `05_test_segmentation_robustness.py`
- `06_generate_clone_atlas.py`
- `07_test_structure_sensitivity_excluding_cycling.R`

## `14_spatial_interactions_and_barriers`

Barrier topology and candidate spatial interactions

Deterministic: `True`

Scripts:

- `00_define_sender_receiver_pairs.py`
- `01_filter_ligand_receptor_database_to_panel.py`
- `02_compute_spatially_constrained_scores.py`
- `03_model_barrier_topology_by_structure.R`
- `04_analyse_barrier_pathways.py`
- `05_prioritise_testable_interactions.py`
- `06_benchmark_against_published_barrier_studies.R`
- `07_ablate_covariates_for_barrier_effect.R`

## `15_hpv_stratified_analysis`

HPV-stratified analysis (consolidated, capped)

Deterministic: `True`

Scripts:

- `00_validate_hpv_metadata_and_probe_signal.py`
- `01_prespecify_primary_hpv_contrasts.py`
- `02_run_prospective_power_simulation.R`
- `03_compare_cellular_composition_patient_level.R`
- `04_compare_ecosystem_and_clone_structure.R`
- `05_run_small_sample_robustness_checks.R`
- `06_prepare_hpv_claim_strength_table.py`

## `16_external_validation_and_generalisation`

Full external validation, triangulation and generalisability

Deterministic: `True`

Scripts:

- `00_define_validation_claims.py`
- `01_acquire_independent_spatial_dataset.py`
- `02_acquire_hnscc_scrna_references.py`
- `03_validate_cell_state_signatures.R`
- `04_validate_ecosystem_signatures_in_bulk.py`
- `05_validate_framework_on_independent_dataset.py`
- `06_compare_with_source_paper_results.py`
- `07_generate_evidence_matrix.py`
- `08_acquire_second_independent_spatial_dataset.py`
- `09_validate_framework_on_second_cancer_type.py`

## `17_statistical_closure_and_release`

Statistical closure, figures, manuscript and reproducible release

Deterministic: `True`

Scripts:

- `00_freeze_primary_results.py`
- `01_control_multiplicity_and_report_effects.R`
- `02_generate_all_main_figures.py`
- `03_generate_all_supplementary_figures.py`
- `04_generate_results_tables.py`
- `05_render_reproducible_manuscript.qmd`
- `05b_render_manuscript_pdf.py`
- `06_build_public_data_release.py`
- `07_build_documented_software_package.py`
- `08_run_end_to_end_reproduction_test.sh`
- `09_run_null_model_calibration_regression.py`
- `11_build_redesigned_manuscript_figures.py`
