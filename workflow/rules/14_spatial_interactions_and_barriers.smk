"""Spatial Interactions and Barriers (14_spatial_interactions_and_barriers): Barrier topology and candidate spatial interactions."""

rule phase14_00_define_sender_receiver_pairs:
    """14.00 -- Predeclares biologically motivated tumour-T cell, fibroblast-T cell, myeloid-T cell and APC-T cell comparisons. Primary output: config/sender_receiver_pairs.yaml."""
    input: "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/00_define_sender_receiver_pairs.done")
    shell: "python3 scripts/14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py"


rule phase14_01_filter_ligand_receptor_database_to_panel:
    """14.01 -- Restricts interaction databases to assayed genes and labels incomplete pairs to prevent overinterpretation. Primary output: references/panel_supported_interactions.tsv."""
    input: "results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/01_filter_ligand_receptor_database_to_panel.done")
    shell: "python3 scripts/14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py"


rule phase14_02_compute_spatially_constrained_scores:
    """14.02 -- Scores candidate-interaction support only across graph-connected cells and compares against degree-preserving nulls. Primary output: data/derived/spatial_interaction_scores.parquet."""
    input:
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/00_define_sender_receiver_pairs.done",
        "results/logs/14_spatial_interactions_and_barriers/.sentinels/01_filter_ligand_receptor_database_to_panel.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/02_compute_spatially_constrained_scores.done")
    shell: "python3 scripts/14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py"


rule phase14_03_model_barrier_topology_by_structure:
    """14.03 -- Tests whether the confirmed clone structure associates with checkpoint, chemokine, TGF-beta, interferon and antigen-presentation programs at barrier interfaces. Primary output: reports/interactions/barrier_topology_models.pdf."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/01_compute_clone_cell_state_composition.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/03_segment_tissue_domains.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/03_model_barrier_topology_by_structure.done")
    shell: "Rscript scripts/14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R"


rule phase14_04_analyse_barrier_pathways:
    """14.04 -- Evaluates programs specifically at spatial interfaces separating excluded clones from tumour. Primary output: data/derived/barrier_pathways.parquet."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/04_analyse_barrier_pathways.done")
    shell: "python3 scripts/14_spatial_interactions_and_barriers/04_analyse_barrier_pathways.py"


rule phase14_05_prioritise_testable_interactions:
    """14.05 -- Ranks candidate interactions by effect size, spatial specificity, cross-patient consistency, panel completeness and external support. Never uses the word 'mechanism'. Primary output: results/interaction_priority_table.tsv."""
    input: "results/logs/14_spatial_interactions_and_barriers/.sentinels/02_compute_spatially_constrained_scores.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/05_prioritise_testable_interactions.done")
    shell: "python3 scripts/14_spatial_interactions_and_barriers/05_prioritise_testable_interactions.py"


rule phase14_06_benchmark_against_published_barrier_studies:
    """14.06 -- Quantitatively compares barrier-topology effect sizes against Grout/Kim 2022 and Hwang et al. 2022. Primary output: reports/interactions/literature_benchmark.pdf."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/06_benchmark_against_published_barrier_studies.done")
    shell: "Rscript scripts/14_spatial_interactions_and_barriers/06_benchmark_against_published_barrier_studies.R"


rule phase14_07_ablate_covariates_for_barrier_effect:
    """14.07 -- Post-hoc follow-up, not part of the original 18-stage plan (see docs/analysis_amendments.md): structured covariate-ablation analysis decomposing which adjustment covariate(s) unmask the significant suppressive-myeloid barrier effect from the weak, non-significant raw bivariate correlation. Primary output: reports/interactions/barrier_covariate_ablation.pdf."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/01_compute_clone_cell_state_composition.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/03_segment_tissue_domains.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done"
    output: touch("results/logs/14_spatial_interactions_and_barriers/.sentinels/07_ablate_covariates_for_barrier_effect.done")
    shell: "Rscript scripts/14_spatial_interactions_and_barriers/07_ablate_covariates_for_barrier_effect.R"
