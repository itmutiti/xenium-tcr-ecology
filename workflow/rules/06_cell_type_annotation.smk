"""Cell Type Annotation (06_cell_type_annotation): Hierarchical cell annotation with uncertainty."""

rule phase06_00_compile_marker_and_reference_registry:
    """06.00 -- Builds versioned marker sets and maps panel genes to feasible cell identities given the targeted assay. Primary output: references/cell_type_marker_registry.tsv."""
    input: "results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/00_compile_marker_and_reference_registry.done")
    shell: "python3 scripts/06_cell_type_annotation/00_compile_marker_and_reference_registry.py"


rule phase06_01_cluster_within_patient_and_jointly:
    """06.01 -- Generates multiple resolutions for exploratory structure while avoiding the assumption that clusters equal cell types. Primary output: data/derived/clustering_assignments.parquet."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/00_separate_gene_and_control_features.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/01_cluster_within_patient_and_jointly.done")
    shell: "python3 scripts/06_cell_type_annotation/01_cluster_within_patient_and_jointly.py"


rule phase06_02_score_major_lineages:
    """06.02 -- Assigns lineage scores for epithelial, T, B, myeloid, fibroblast, endothelial and other compartments. Primary output: data/derived/lineage_scores.parquet."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/00_compile_marker_and_reference_registry.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/02_score_major_lineages.done")
    shell: "python3 scripts/06_cell_type_annotation/02_score_major_lineages.py"


rule phase06_03_map_external_scrna_reference:
    """06.03 -- Maps a suitable HNSCC single-cell reference with label transfer restricted to overlapping genes, and reports the transfer-confidence degradation caused by the restricted gene overlap. Primary output: data/derived/reference_labels.parquet."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/08_acquire_companion_scrna_and_vdj_reference.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/03_map_external_scrna_reference.done")
    shell: "python3 scripts/06_cell_type_annotation/03_map_external_scrna_reference.py"


rule phase06_04_resolve_t_cell_substates:
    """06.04 -- Annotates CD4, CD8, Treg, cycling, cytotoxic, exhausted/dysfunctional and ambiguous T-cell states. Primary output: data/derived/t_cell_states.parquet."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/03_calculate_program_scores.done",
        "results/logs/06_cell_type_annotation/.sentinels/02_score_major_lineages.done",
        "results/logs/06_cell_type_annotation/.sentinels/03_map_external_scrna_reference.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done")
    shell: "Rscript scripts/06_cell_type_annotation/04_resolve_t_cell_substates.R"


rule phase06_05_resolve_myeloid_and_stromal_substates:
    """06.05 -- Defines macrophage, DC, fibroblast, endothelial and perivascular subtypes at the resolution the panel supports. Primary output: data/derived/tme_substates.parquet."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/02_score_major_lineages.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/05_resolve_myeloid_and_stromal_substates.done")
    shell: "Rscript scripts/06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R"


rule phase06_06_integrate_annotation_evidence:
    """06.06 -- Combines marker, reference, cluster and spatial evidence into labels plus calibrated confidence and ambiguity flags. Primary output: data/derived/final_cell_annotations.parquet."""
    input:
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/01_cluster_within_patient_and_jointly.done",
        "results/logs/06_cell_type_annotation/.sentinels/02_score_major_lineages.done",
        "results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done",
        "results/logs/06_cell_type_annotation/.sentinels/05_resolve_myeloid_and_stromal_substates.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done")
    shell: "python3 scripts/06_cell_type_annotation/06_integrate_annotation_evidence.py"


rule phase06_07_blinded_annotation_review:
    """06.07 -- Creates blinded spatial panels for expert review and records adjudications without overwriting original predictions. Primary output: reports/annotation/adjudication_log.tsv."""
    input:
        "results/logs/.checkpoints/canonical_objects.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/06_cell_type_annotation/.sentinels/07_blinded_annotation_review.done")
    shell: "python3 scripts/06_cell_type_annotation/07_blinded_annotation_review.py"


rule phase06_08_acquire_companion_scrna_and_vdj_reference:
    """06.08 -- Downloads and verifies GSE287301 (McCord et al. 2026's own companion scRNA-seq + paired VDJ dataset): public, unauthenticated NCBI GEO download, no manual step required. Primary output: data/external/GSE287301/.
    Real, direct dependents (not numeric proximity): `phase06_03_map_
    external_scrna_reference` (this rule's own `input:`) and
    `08_tcr_clonal_analysis/phase08_09_validate_probe_clones_against_
    paired_vdj_ground_truth` (that rule's own `input:`, in a different
    phase). Numbered as a Cell Type Annotation milestone -- its first
    consumer -- rather than renumbered into an earlier position, to avoid
    disturbing any existing rule's identity."""
    output: touch("results/logs/06_cell_type_annotation/.sentinels/08_acquire_companion_scrna_and_vdj_reference.done")
    shell: "python3 scripts/06_cell_type_annotation/08_acquire_companion_scrna_and_vdj_reference.py"
