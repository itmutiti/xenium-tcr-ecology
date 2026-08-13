"""Niche and Ecosystem Discovery (10_niche_and_ecosystem_discovery): Discovery of spatial ecosystems and tissue domains."""

rule phase10_00_compute_cell_type_neighbourhood_enrichment:
    """10.00 -- Tests pairwise cell-type adjacency using within-section constrained permutations (calibrated in Spatial Graph Construction and Calibration). Primary output: data/derived/neighbourhood_enrichment.parquet."""
    input: "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/00_compute_cell_type_neighbourhood_enrichment.done")
    shell: "python3 scripts/10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py"


rule phase10_01_compute_local_neighbourhood_compositions:
    """10.01 -- Creates cell-centred neighbourhood composition vectors at several spatial scales. Primary output: data/derived/local_compositions.parquet."""
    input:
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/01_prune_graphs_for_tissue_gaps.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done")
    shell: "python3 scripts/10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py"


rule phase10_02_discover_neighbourhood_archetypes:
    """10.02 -- Identifies recurrent neighbourhoods with consensus clustering or topic modelling and quantifies stability. Primary output: data/derived/neighbourhood_archetypes.parquet."""
    input: "results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/02_discover_neighbourhood_archetypes.done")
    shell: "Rscript scripts/10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R"


rule phase10_03_segment_tissue_domains:
    """10.03 -- Forms contiguous domains from neighbourhood states while preserving boundaries and rare niches. Primary output: data/derived/tissue_domains.parquet."""
    input:
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/02_discover_neighbourhood_archetypes.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/03_segment_tissue_domains.done")
    shell: "python3 scripts/10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py"


rule phase10_04_annotate_ecosystems_with_blinded_rules:
    """10.04 -- Assigns descriptive ecosystem labels only after unsupervised discovery, using a documented rubric. Primary output: metadata/ecosystem_annotation.tsv."""
    input:
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/02_discover_neighbourhood_archetypes.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/04_annotate_ecosystems_with_blinded_rules.done")
    shell: "python3 scripts/10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py"


rule phase10_05_quantify_ecosystem_abundance_and_topology:
    """10.05 -- Measures abundance, fragmentation, interface length, mixing, compactness and relation to tumour borders. Primary output: data/derived/ecosystem_metrics.parquet."""
    input:
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/03_segment_tissue_domains.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/04_annotate_ecosystems_with_blinded_rules.done",
        "results/logs/09_spatial_graph_construction_and_calibration/.sentinels/03_construct_primary_cell_graph.done",
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/05_extract_tumour_boundaries.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/05_quantify_ecosystem_abundance_and_topology.done")
    shell: "python3 scripts/10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py"


rule phase10_06_test_patient_recurrence:
    """10.06 -- Fits patient-level models with technical replicates nested within patient; HPV is not tested here. Primary output: reports/niches/recurrence_models.pdf."""
    input: "results/logs/10_niche_and_ecosystem_discovery/.sentinels/05_quantify_ecosystem_abundance_and_topology.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/06_test_patient_recurrence.done")
    shell: "Rscript scripts/10_niche_and_ecosystem_discovery/06_test_patient_recurrence.R"


rule phase10_07_leave_one_patient_out_niche_stability:
    """10.07 -- Tests whether archetypes and labels remain identifiable when each patient is withheld. Primary output: reports/niches/LOPO_stability.pdf."""
    input:
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/01_compute_local_neighbourhood_compositions.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/02_discover_neighbourhood_archetypes.done"
    output: touch("results/logs/10_niche_and_ecosystem_discovery/.sentinels/07_leave_one_patient_out_niche_stability.done")
    shell: "python3 scripts/10_niche_and_ecosystem_discovery/07_leave_one_patient_out_niche_stability.py"
