"""Tumour Epithelium Characterisation (07_tumour_epithelium_characterisation): Malignant epithelial-state inference and tumour-boundary reconstruction."""

rule phase07_00_subset_and_recluster_epithelial_cells:
    """07.00 -- Characterises epithelial heterogeneity separately from immune/stromal structure. Primary output: data/objects/epithelial_subset.h5ad."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/00_subset_and_recluster_epithelial_cells.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/00_subset_and_recluster_epithelial_cells.py"


rule phase07_01_score_malignancy_and_normal_epithelium:
    """07.01 -- Combines tumour keratins, HPV genes, stress/EMT programs and reference mapping to estimate malignancy probability. Primary output: data/derived/malignancy_scores.parquet."""
    input: "results/logs/07_tumour_epithelium_characterisation/.sentinels/00_subset_and_recluster_epithelial_cells.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/01_score_malignancy_and_normal_epithelium.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/01_score_malignancy_and_normal_epithelium.py"


rule phase07_02_cross_validate_against_morphology:
    """07.02 -- Cross-checks the transcriptional malignancy call against any co-registered morphology/pathology annotation and reports concordance explicitly, to bound circularity risk. Primary output: reports/tumour/morphology_concordance.pdf."""
    input:
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/01_score_malignancy_and_normal_epithelium.done",
        "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/02_cross_validate_against_morphology.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/02_cross_validate_against_morphology.py"


rule phase07_03_infer_cnv_appendix_only:
    """07.03 -- Attempts CNV-like inference from the 477-feature panel; restricted to an appendix sensitivity output and barred from any primary claim path. Primary output: reports/tumour/cnv_appendix.pdf."""
    input: "results/logs/07_tumour_epithelium_characterisation/.sentinels/01_score_malignancy_and_normal_epithelium.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/03_infer_cnv_appendix_only.done")
    shell: "Rscript scripts/07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R"


rule phase07_04_construct_tumour_region_masks:
    """07.04 -- Forms spatially coherent tumour regions from malignant-cell probabilities and removes isolated false positives. Primary output: data/derived/tumour_masks/."""
    input: "results/logs/07_tumour_epithelium_characterisation/.sentinels/01_score_malignancy_and_normal_epithelium.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/04_construct_tumour_region_masks.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py"


rule phase07_05_extract_tumour_boundaries:
    """07.05 -- Creates tumour-border geometries and signed distance-to-boundary values for all cells. Primary output: data/derived/tumour_boundaries.parquet."""
    input: "results/logs/07_tumour_epithelium_characterisation/.sentinels/04_construct_tumour_region_masks.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/05_extract_tumour_boundaries.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py"


rule phase07_06_validate_boundaries_against_morphology:
    """07.06 -- Overlays inferred boundaries on morphology images and quantifies manual-review agreement. Primary output: reports/tumour/boundary_validation.pdf."""
    input:
        "results/logs/07_tumour_epithelium_characterisation/.sentinels/04_construct_tumour_region_masks.done",
        "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/06_validate_boundaries_against_morphology.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/06_validate_boundaries_against_morphology.py"


rule phase07_07_define_invasive_front_and_compartments:
    """07.07 -- Labels tumour core, inner margin, outer margin and distal stroma using predeclared distance bands and sensitivity analyses. Primary output: data/derived/spatial_compartments.parquet."""
    input: "results/logs/07_tumour_epithelium_characterisation/.sentinels/05_extract_tumour_boundaries.done"
    output: touch("results/logs/07_tumour_epithelium_characterisation/.sentinels/07_define_invasive_front_and_compartments.done")
    shell: "python3 scripts/07_tumour_epithelium_characterisation/07_define_invasive_front_and_compartments.py"
