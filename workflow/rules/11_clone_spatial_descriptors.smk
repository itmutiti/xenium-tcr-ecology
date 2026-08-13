"""Clone Spatial Descriptors (11_clone_spatial_descriptors): Clone-level ecological descriptors and structure test (provisional)."""

rule phase11_00_compute_clone_spatial_descriptors_rarefied:
    """11.00 -- Calculates dispersion, clustering, tumour distance, border enrichment and domain occupancy using rarefaction/size-normalisation. Primary output: data/derived/clone_spatial_descriptors.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/05_extract_tumour_boundaries.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/07_define_invasive_front_and_compartments.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/03_segment_tissue_domains.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/00_compute_clone_spatial_descriptors_rarefied.done")
    shell: "python3 scripts/11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py"


rule phase11_01_compute_clone_cell_state_composition:
    """11.01 -- Summarises cytotoxic, exhausted, activated, proliferative and ambiguous states per clone with uncertainty. Primary output: data/derived/clone_state_composition.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/01_compute_clone_cell_state_composition.done")
    shell: "python3 scripts/11_clone_spatial_descriptors/01_compute_clone_cell_state_composition.py"


rule phase11_02_quantify_clone_tumour_engagement:
    """11.02 -- Measures malignant-cell adjacency, contact normalised by opportunity, penetration and interface localisation; cross-checked against Quality Control spillover-risk and resegmentation outputs. Primary output: data/derived/clone_tumour_engagement.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/05_extract_tumour_boundaries.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/00_compute_clone_spatial_descriptors_rarefied.done",
        "results/logs/04_quality_control/.sentinels/04_estimate_transcript_spillover.done",
        "results/logs/04_quality_control/.sentinels/05_resegment_reference_subset.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done")
    shell: "python3 scripts/11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py"


rule phase11_03_quantify_clone_apc_support:
    """11.03 -- Measures proximity and neighbourhood enrichment with dendritic cells, macrophages and antigen-presentation programs. Primary output: data/derived/clone_apc_support.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/03_quantify_clone_apc_support.done")
    shell: "python3 scripts/11_clone_spatial_descriptors/03_quantify_clone_apc_support.py"


rule phase11_04_quantify_stromal_and_myeloid_barriers:
    """11.04 -- Builds clone-specific fibroblast, vascular and suppressive-myeloid barrier metrics along paths to tumour, using calibrated null models. Primary output: data/derived/clone_barrier_metrics.parquet."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done")
    shell: "python3 scripts/11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py"


rule phase11_05_test_discrete_vs_continuous_structure:
    """11.05 -- Compares a discrete clustering model against a continuous archetypal/latent-trajectory model; mandatory before any taxonomy step. Primary output: reports/clone_ecology/structure_test.pdf."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/00_compute_clone_spatial_descriptors_rarefied.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/01_compute_clone_cell_state_composition.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/03_quantify_clone_apc_support.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/05_test_discrete_vs_continuous_structure.done")
    shell: "Rscript scripts/11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R"


rule phase11_06_discover_provisional_structure:
    """11.06 -- Fits whichever structure 11.05 supports (discrete fates or a continuous score), without imposing categories the data don't show. Primary output: data/derived/clone_ecological_structure.parquet."""
    input: "results/logs/11_clone_spatial_descriptors/.sentinels/05_test_discrete_vs_continuous_structure.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/06_discover_provisional_structure.done")
    shell: "Rscript scripts/11_clone_spatial_descriptors/06_discover_provisional_structure.R"


rule phase11_07_freeze_provisional_taxonomy_version:
    """11.07 -- Freezes taxonomy/structure definitions as taxonomy_version=v1 and routes to External Checkpoint Validation for an external sanity check. Primary output: governance/taxonomy_version_log.tsv."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/05_test_discrete_vs_continuous_structure.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/06_discover_provisional_structure.done"
    output: touch("results/logs/11_clone_spatial_descriptors/.sentinels/07_freeze_provisional_taxonomy_version.done")
    shell: "python3 scripts/11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py"
