"""Spatial Graph Construction and Calibration (09_spatial_graph_construction_and_calibration): Spatial graph construction, sensitivity framework and null-model calibration."""

rule phase09_00_generate_candidate_spatial_graphs:
    """09.00 -- Builds radius, k-nearest-neighbour, Delaunay and boundary-contact graphs for every section. Primary output: data/graphs/candidate/."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/00_generate_candidate_spatial_graphs.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py"


rule phase09_01_prune_graphs_for_tissue_gaps:
    """09.01 -- Removes Delaunay/radius edges that bridge tissue folds, holes or QC-excluded regions. Primary output: data/graphs/pruned/."""
    input: "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/00_generate_candidate_spatial_graphs.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/01_prune_graphs_for_tissue_gaps.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py"


rule phase09_02_calibrate_graph_parameters:
    """09.02 -- Selects biologically interpretable scales using cell size, interaction ranges, graph connectivity and stability. Primary output: config/graph_parameters.yaml."""
    input: "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/01_prune_graphs_for_tissue_gaps.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/02_calibrate_graph_parameters.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py"


rule phase09_03_construct_primary_cell_graph:
    """09.03 -- Creates the prespecified primary graph with node/edge metadata and patient-separated components. Primary output: data/graphs/primary_graphs/."""
    input:
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/02_calibrate_graph_parameters.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py"


rule phase09_04_construct_tumour_tcell_bipartite_graph:
    """09.04 -- Represents direct and near-direct malignant-cell/T-cell relationships for contact-focused analyses. Primary output: data/graphs/tumour_tcell/."""
    input:
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/01_prune_graphs_for_tissue_gaps.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/04_construct_tumour_region_masks.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/04_construct_tumour_tcell_bipartite_graph.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/04_construct_tumour_tcell_bipartite_graph.py"


rule phase09_05_construct_clone_induced_subgraphs:
    """09.05 -- Extracts all cells belonging to each clone plus local microenvironment shells. Primary output: data/graphs/clones/."""
    input:
        "results/logs/08_tcr_clonal_analysis/.sentinels/08_generate_tcr_release_report.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/05_construct_clone_induced_subgraphs.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/05_construct_clone_induced_subgraphs.py"


rule phase09_06_run_graph_sensitivity_grid:
    """09.06 -- Repeats core metrics over plausible radii/k values to identify robust conclusions. Primary output: reports/graphs/sensitivity_grid.parquet."""
    input:
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/01_prune_graphs_for_tissue_gaps.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/04_construct_tumour_region_masks.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/06_run_graph_sensitivity_grid.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/06_run_graph_sensitivity_grid.py"


rule phase09_07_generate_synthetic_ground_truth_patterns:
    """09.07 -- Simulates spatial point patterns with known, controllable effect sizes and known clone-spatial relationships. Primary output: data/graphs/synthetic/."""
    input: "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/02_calibrate_graph_parameters.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/07_generate_synthetic_ground_truth_patterns.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py"


rule phase09_08_run_calibration_suite_on_synthetic_data:
    """09.08 -- Confirms each planned null model achieves nominal Type I error and adequate power at n=10-patient scale before it is applied to observed data. Primary output: reports/graphs/null_model_calibration.pdf."""
    input: "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/07_generate_synthetic_ground_truth_patterns.done"
    output: touch("results/logs/09_spatial_graph_construction_and_calibration/.sentinels/08_run_calibration_suite_on_synthetic_data.done")
    shell: "python3 scripts/09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py"
