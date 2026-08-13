"""Clone Ecology Confirmatory Models (13_clone_ecology_confirmatory_models): Confirmatory variance-partitioning and clone ecological modelling."""

rule phase13_00_load_frozen_taxonomy_version:
    """13.00 -- Loads the exact frozen taxonomy_version from External Checkpoint Validation; refuses to run if any upstream input hash has changed. Primary output: data/releases/v1_clone_structure/."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/07_freeze_provisional_taxonomy_version.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/03_decide_freeze_or_revise.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done")
    shell: "python3 scripts/13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py"


rule phase13_01_fit_variance_partition_models:
    """13.01 -- Fits hierarchical/Bayesian models decomposing variance in clone spatial/phenotypic descriptors into identity, context and patient components. Primary output: reports/clone_ecology/variance_partition_models.pdf."""
    input: "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/01_fit_variance_partition_models.done")
    shell: "Rscript scripts/13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R"


rule phase13_02_fit_hierarchical_clone_models:
    """13.02 -- Fits mixed or Bayesian models for any confirmed categorical structure, nested in patient and section. Primary output: reports/clone_ecology/hierarchical_models.pdf."""
    input: "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/02_fit_hierarchical_clone_models.done")
    shell: "Rscript scripts/13_clone_ecology_confirmatory_models/02_fit_hierarchical_clone_models.R"


rule phase13_03_test_clone_size_as_confounder:
    """13.03 -- Determines which associations survive after matching or adjustment for clone size and detection power. Primary output: reports/clone_ecology/clone_size_sensitivity.pdf."""
    input: "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/03_test_clone_size_as_confounder.done")
    shell: "Rscript scripts/13_clone_ecology_confirmatory_models/03_test_clone_size_as_confounder.R"


rule phase13_04_run_leave_one_patient_out_stability:
    """13.04 -- Confirms the variance partition and any categorical structure are not driven by a single patient. Primary output: reports/clone_ecology/lopo_stability.pdf."""
    input: "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/04_run_leave_one_patient_out_stability.done")
    shell: "Rscript scripts/13_clone_ecology_confirmatory_models/04_run_leave_one_patient_out_stability.R"


rule phase13_05_test_segmentation_robustness:
    """13.05 -- Re-runs the confirmatory models on the Quality Control resegmented reference subset and reports whether conclusions hold. Primary output: reports/clone_ecology/segmentation_robustness.pdf."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done",
        "results/logs/04_quality_control/.sentinels/05_resegment_reference_subset.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/05_test_segmentation_robustness.done")
    shell: "python3 scripts/13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py"


rule phase13_06_generate_clone_atlas:
    """13.06 -- Produces standardised maps and summaries for every high-confidence clone, preventing example-selection bias. Primary output: reports/clone_ecology/clone_atlas.html."""
    input:
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done",
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/01_compute_clone_cell_state_composition.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/02_quantify_clone_tumour_engagement.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/03_quantify_clone_apc_support.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/04_quantify_stromal_and_myeloid_barriers.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/06_generate_clone_atlas.done")
    shell: "python3 scripts/13_clone_ecology_confirmatory_models/06_generate_clone_atlas.py"


rule phase13_07_test_structure_sensitivity_excluding_cycling:
    """13.07 -- Post-hoc sensitivity check, not part of the original 18-stage plan (see docs/analysis_amendments.md): re-derives the continuous ecological-structure score and variance partition with cycling_fraction excluded, to test how much Q2's headline result depends on the one externally-flagged input feature. Does not touch the frozen v1_provisional release. Primary output: reports/clone_ecology/structure_sensitivity_excluding_cycling.pdf."""
    input:
        "results/logs/11_clone_spatial_descriptors/.sentinels/06_discover_provisional_structure.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/01_fit_variance_partition_models.done"
    output: touch("results/logs/13_clone_ecology_confirmatory_models/.sentinels/07_test_structure_sensitivity_excluding_cycling.done")
    shell: "Rscript scripts/13_clone_ecology_confirmatory_models/07_test_structure_sensitivity_excluding_cycling.R"
