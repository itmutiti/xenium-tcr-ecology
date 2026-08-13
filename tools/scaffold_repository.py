#!/usr/bin/env python3
"""
Idempotent repository scaffold generator for hnscc-xenium-tcr-ecology.

Creates the full directory tree, phase script stubs, workflow rule
skeletons, and a governance/policy verification pass, per the project's
18-phase specification (src/ package split, phase-first script folders,
governance-as-code checks). Re-running it is safe: existing files are
left untouched unless --force is passed, and it ends with a policy
verification pass that fails loudly on violation.

`PHASES`/`CHECKPOINTS` below are kept in sync with the `scripts/`/
`workflow/rules/` folder names (indices 1-17); treat `workflow/rules/*.smk`
and `scripts/` as the authoritative current structure and this file as
the generator for extending it consistently -- e.g. the checkpoint
chain's final rule is `computational_analysis_complete`, with a
separate, explicitly-named `publication_release` target beyond it.

Usage:
    python3 tools/scaffold_repository.py [--project-root PATH] [--force]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Phase / script specification.
# (phase_index, folder_name, purpose) and, per phase, an ordered list of
# (order_code, filename, language, what_it_accomplishes, primary_output).
# ---------------------------------------------------------------------------

PHASES = [
    (
        1,
        "01_project_setup_and_governance",
        "Research governance, computational environment and reproducibility",
        [
            (
                "01.00",
                "00_validate_project_scope.py",
                "py",
                "Creates a machine-readable project charter recording the primary research questions; separates published findings from original hypotheses.",
                "project_charter.yaml",
            ),
            (
                "01.01",
                "01_build_sample_manifest.py",
                "py",
                "Parses GEO sample metadata into one canonical manifest with patient, run, replicate, HPV status, specimen type, file paths and checksums, and records the exact GEO record revision date used.",
                "metadata/sample_manifest.tsv",
            ),
            (
                "01.02",
                "02_generate_data_dictionary.py",
                "py",
                "Defines every metadata field, unit, allowed value, missingness code and derivation rule.",
                "metadata/data_dictionary.xlsx",
            ),
            (
                "01.03",
                "03_lock_software_environments.sh",
                "sh",
                "Builds pinned conda environments for Python and R; records OS and compiler information.",
                "env/environment.yml; env/renv.lock",
            ),
            (
                "01.04",
                "04_build_container_images.sh",
                "sh",
                "Builds a Docker/Apptainer image pinned to the locked environments, for the eventual frozen release.",
                "containers/release.sif; containers/Dockerfile",
            ),
            (
                "01.05",
                "05_initialise_reproducible_workflow.py",
                "py",
                "Initialises Snakemake profiles, configuration templates, log directories, benchmark reporting and deterministic seeds.",
                "Snakefile; config/config.yaml; profiles/",
            ),
            (
                "01.06",
                "06_create_analysis_registry.py",
                "py",
                "Registers each planned analysis with hypothesis, unit of analysis, exclusion criteria, primary endpoint, multiplicity family and taxonomy_version tracking.",
                "governance/analysis_registry.tsv",
            ),
        ],
    ),
    (
        2,
        "02_raw_data_ingestion",
        "Data acquisition, integrity verification and immutable staging",
        [
            (
                "02.00",
                "00_query_geo_accession.py",
                "py",
                "Queries GSE300147 metadata and records accession status, sample list, file sizes, retrieval date and record revision.",
                "metadata/geo_snapshot.json",
            ),
            (
                "02.01",
                "01_download_geo_raw_archive.sh",
                "sh",
                "Downloads the raw GEO archive with resumable transfer and conservative retry logic.",
                "data/raw/GSE300147_RAW.tar",
            ),
            (
                "02.02",
                "02_verify_archive_checksums.py",
                "py",
                "Calculates SHA-256 checksums, verifies file completeness, produces an immutable manifest.",
                "data/raw/SHA256SUMS; reports/integrity.tsv",
            ),
            (
                "02.03",
                "03_extract_archive_safely.py",
                "py",
                "Extracts without path traversal, preserves timestamps, maps files to GEO sample identifiers.",
                "data/staged/<sample_id>/",
            ),
            (
                "02.04",
                "04_inventory_xenium_files.py",
                "py",
                "Detects H5, Parquet, TIFF and auxiliary files; flags missing mandatory assets and format-version differences.",
                "metadata/xenium_file_inventory.tsv",
            ),
            (
                "02.05",
                "05_standardise_sample_directory_layout.py",
                "py",
                "Creates a consistent per-section directory structure without modifying source files, using links where possible.",
                "data/standardised/<section_id>/",
            ),
            (
                "02.06",
                "06_freeze_raw_data_snapshot.sh",
                "sh",
                "Marks raw and staged data read-only and records provenance and storage usage.",
                "reports/raw_data_freeze.txt",
            ),
        ],
    ),
    (
        3,
        "03_spatialdata_import",
        "Xenium import, coordinate harmonisation and canonical data objects",
        [
            (
                "03.00",
                "00_detect_xenium_format_version.py",
                "py",
                "Inspects file schemas, chooses compatible readers, cross-validates against a second independent reader implementation, records deviations.",
                "metadata/format_versions.tsv",
            ),
            (
                "03.01",
                "01_import_each_section_to_spatialdata.py",
                "py",
                "Loads counts, cell metadata, transcript coordinates, boundaries and morphology into one SpatialData object per section.",
                "data/objects/spatialdata/<section>.zarr",
            ),
            (
                "03.02",
                "02_validate_coordinate_systems.py",
                "py",
                "Confirms micrometre units, image orientation, origin and alignment among transcripts, centroids, polygons and images; confirms polygon-vs-centroid convention consistently across sections.",
                "reports/coordinate_validation/",
            ),
            (
                "03.03",
                "03_create_anndata_expression_objects.py",
                "py",
                "Extracts a clean cell-by-gene AnnData object while preserving links to spatial elements.",
                "data/objects/anndata/<section>.h5ad",
            ),
            (
                "03.04",
                "04_attach_clinical_and_technical_metadata.py",
                "py",
                "Joins patient and run metadata using validated keys and prevents silent many-to-many joins.",
                "updated h5ad/zarr objects; join audit",
            ),
            (
                "03.05",
                "05_build_combined_analysis_object.py",
                "py",
                "Concatenates sections with globally unique cell identifiers and explicit patient/run hierarchy.",
                "data/objects/hnscc_xenium_combined.h5ad",
            ),
            (
                "03.06",
                "06_export_r_interoperability_objects.py",
                "py",
                "Creates Seurat-compatible and parquet exports for R analyses without losing spatial coordinates.",
                "data/objects/r_exports/",
            ),
        ],
    ),
    (
        4,
        "04_quality_control",
        "Multilevel quality control, principled exclusion and transcript-integrity robustness",
        [
            (
                "04.00",
                "00_compute_cell_level_qc_metrics.py",
                "py",
                "Calculates transcript counts, detected genes, cell/nucleus area, control ratios and transcript density.",
                "data/derived/cell_qc_metrics.parquet",
            ),
            (
                "04.01",
                "01_compute_transcript_level_qc_metrics.py",
                "py",
                "Profiles Q-values, unassigned transcripts, controls, nuclear overlap and spatial transcript-density artefacts.",
                "data/derived/transcript_qc_metrics.parquet",
            ),
            (
                "04.02",
                "02_detect_spatial_qc_artifacts.py",
                "py",
                "Finds damaged regions, edge effects, holes, striping and local decoding failures using spatial statistics.",
                "reports/qc/spatial_artifact_masks/",
            ),
            (
                "04.03",
                "03_assess_segmentation_quality.py",
                "py",
                "Measures cell geometry, nucleus containment, transcript assignment and suspicious multinucleation; samples cells for visual review.",
                "reports/qc/segmentation_review.html",
            ),
            (
                "04.04",
                "04_estimate_transcript_spillover.py",
                "py",
                "Flags cells near segmentation boundaries adjacent to a different predicted cell type; estimates a per-cell spillover-risk score using distance-to-boundary and neighbour-identity weighting.",
                "data/derived/spillover_risk.parquet",
            ),
            (
                "04.05",
                "05_resegment_reference_subset.py",
                "py",
                "Independently resegments a representative subset of sections with an alternative transcript-reassignment method; used solely to test whether headline spatial results are robust to segmentation choice.",
                "data/objects/resegmented_subset/",
            ),
            (
                "04.06",
                "06_define_qc_thresholds_hierarchically.R",
                "R",
                "Uses robust, section-aware thresholds rather than global cut-offs; documents sensitivity alternatives.",
                "config/qc_thresholds.yaml",
            ),
            (
                "04.07",
                "07_apply_qc_filters_with_audit_trail.py",
                "py",
                "Applies filters without deleting data and records a reason code for every excluded cell and transcript.",
                "data/objects/qc_filtered.h5ad; exclusion_log.tsv",
            ),
            (
                "04.08",
                "08_assess_replicate_concordance.R",
                "R",
                "Quantifies technical replicate agreement in counts, cell composition, gene expression and spatial statistics, correctly modelled as day nested within patient.",
                "reports/qc/replicate_concordance.pdf",
            ),
            (
                "04.09",
                "09_generate_qc_release_report.py",
                "py",
                "Produces a section-by-section QC report and a go/no-go decision before biological analysis.",
                "reports/qc/QC_release_report.html",
            ),
        ],
    ),
    (
        5,
        "05_preprocessing_and_normalisation",
        "Expression normalisation, feature engineering and technical variation",
        [
            (
                "05.00",
                "00_separate_gene_and_control_features.py",
                "py",
                "Separates biological genes, negative controls, HPV probes and TCR/CDR3 probes into explicit feature classes.",
                "metadata/feature_annotation.tsv",
            ),
            (
                "05.01",
                "01_construct_analysis_count_layers.py",
                "py",
                "Preserves raw counts and creates normalised, variance-stabilised and binary-detection layers where appropriate.",
                "layers in h5ad",
            ),
            (
                "05.02",
                "02_evaluate_normalisation_strategies.R",
                "R",
                "Compares log-normalisation, Pearson residuals and detection-based approaches using replicate stability and negative-control probe behaviour.",
                "reports/preprocess/normalisation_benchmark.pdf",
            ),
            (
                "05.03",
                "03_calculate_program_scores.py",
                "py",
                "Computes curated cytotoxicity, exhaustion, activation, interferon, proliferation, stress, EMT and antigen-presentation scores.",
                "data/derived/program_scores.parquet",
            ),
            (
                "05.04",
                "04_model_technical_covariates.R",
                "R",
                "Estimates contributions of run, section, depth and control burden without erasing patient biology; run is modelled explicitly as nested within patient.",
                "reports/preprocess/variance_partition.pdf",
            ),
            (
                "05.05",
                "05_create_primary_analysis_matrix.py",
                "py",
                "Freezes the matrix and covariates used by downstream analyses and records its exact hash.",
                "data/releases/v1_primary_analysis/",
            ),
        ],
    ),
    (
        6,
        "06_cell_type_annotation",
        "Hierarchical cell annotation with uncertainty",
        [
            (
                "06.00",
                "00_compile_marker_and_reference_registry.py",
                "py",
                "Builds versioned marker sets and maps panel genes to feasible cell identities given the targeted assay.",
                "references/cell_type_marker_registry.tsv",
            ),
            (
                "06.01",
                "01_cluster_within_patient_and_jointly.py",
                "py",
                "Generates multiple resolutions for exploratory structure while avoiding the assumption that clusters equal cell types.",
                "data/derived/clustering_assignments.parquet",
            ),
            (
                "06.02",
                "02_score_major_lineages.py",
                "py",
                "Assigns lineage scores for epithelial, T, B, myeloid, fibroblast, endothelial and other compartments.",
                "data/derived/lineage_scores.parquet",
            ),
            (
                "06.03",
                "03_map_external_scrna_reference.py",
                "py",
                "Maps a suitable HNSCC single-cell reference with label transfer restricted to overlapping genes, and reports the transfer-confidence degradation caused by the restricted gene overlap.",
                "data/derived/reference_labels.parquet",
            ),
            (
                "06.04",
                "04_resolve_t_cell_substates.R",
                "R",
                "Annotates CD4, CD8, Treg, cycling, cytotoxic, exhausted/dysfunctional and ambiguous T-cell states.",
                "data/derived/t_cell_states.parquet",
            ),
            (
                "06.05",
                "05_resolve_myeloid_and_stromal_substates.R",
                "R",
                "Defines macrophage, DC, fibroblast, endothelial and perivascular subtypes at the resolution the panel supports.",
                "data/derived/tme_substates.parquet",
            ),
            (
                "06.06",
                "06_integrate_annotation_evidence.py",
                "py",
                "Combines marker, reference, cluster and spatial evidence into labels plus calibrated confidence and ambiguity flags.",
                "data/derived/final_cell_annotations.parquet",
            ),
            (
                "06.07",
                "07_blinded_annotation_review.py",
                "py",
                "Creates blinded spatial panels for expert review and records adjudications without overwriting original predictions.",
                "reports/annotation/adjudication_log.tsv",
            ),
        ],
    ),
    (
        7,
        "07_tumour_epithelium_characterisation",
        "Malignant epithelial-state inference and tumour-boundary reconstruction",
        [
            (
                "07.00",
                "00_subset_and_recluster_epithelial_cells.py",
                "py",
                "Characterises epithelial heterogeneity separately from immune/stromal structure.",
                "data/objects/epithelial_subset.h5ad",
            ),
            (
                "07.01",
                "01_score_malignancy_and_normal_epithelium.py",
                "py",
                "Combines tumour keratins, HPV genes, stress/EMT programs and reference mapping to estimate malignancy probability.",
                "data/derived/malignancy_scores.parquet",
            ),
            (
                "07.02",
                "02_cross_validate_against_morphology.py",
                "py",
                "Cross-checks the transcriptional malignancy call against any co-registered morphology/pathology annotation and reports concordance explicitly, to bound circularity risk.",
                "reports/tumour/morphology_concordance.pdf",
            ),
            (
                "07.03",
                "03_infer_cnv_appendix_only.R",
                "R",
                "Attempts CNV-like inference from the 477-feature panel; restricted to an appendix sensitivity output and barred from any primary claim path.",
                "reports/tumour/cnv_appendix.pdf",
            ),
            (
                "07.04",
                "04_construct_tumour_region_masks.py",
                "py",
                "Forms spatially coherent tumour regions from malignant-cell probabilities and removes isolated false positives.",
                "data/derived/tumour_masks/",
            ),
            (
                "07.05",
                "05_extract_tumour_boundaries.py",
                "py",
                "Creates tumour-border geometries and signed distance-to-boundary values for all cells.",
                "data/derived/tumour_boundaries.parquet",
            ),
            (
                "07.06",
                "06_validate_boundaries_against_morphology.py",
                "py",
                "Overlays inferred boundaries on morphology images and quantifies manual-review agreement.",
                "reports/tumour/boundary_validation.pdf",
            ),
            (
                "07.07",
                "07_define_invasive_front_and_compartments.py",
                "py",
                "Labels tumour core, inner margin, outer margin and distal stroma using predeclared distance bands and sensitivity analyses.",
                "data/derived/spatial_compartments.parquet",
            ),
        ],
    ),
    (
        8,
        "08_tcr_clonal_analysis",
        "TCR probe decoding, clone assignment, validation and ascertainment audit",
        [
            (
                "08.00",
                "00_identify_tcr_cdr3_probe_features.py",
                "py",
                "Separates patient-specific CDR3 probes from conventional T-cell genes using feature metadata and naming rules.",
                "metadata/tcr_probe_registry.tsv",
            ),
            (
                "08.01",
                "01_map_tcr_probes_to_patients.py",
                "py",
                "Ensures patient-specific probes are only evaluated in their intended specimens and detects leakage or naming conflicts.",
                "reports/tcr/patient_probe_audit.tsv",
            ),
            (
                "08.02",
                "02_document_clone_ascertainment.py",
                "py",
                "Records how each probed clonotype was selected relative to each patient's full repertoire and publishes this as an explicit boundary condition on every downstream generalisability claim.",
                "metadata/clone_ascertainment.tsv",
            ),
            (
                "08.03",
                "03_call_cell_level_tcr_detections.py",
                "py",
                "Calls clone detections from transcript counts with explicit thresholds, controls and multi-probe ambiguity handling.",
                "data/derived/tcr_cell_calls.parquet",
            ),
            (
                "08.04",
                "04_estimate_false_positive_tcr_calls.R",
                "R",
                "Uses off-patient probes, non-T cells and spatial permutation as empirical negative controls.",
                "reports/tcr/false_positive_model.pdf",
            ),
            (
                "08.05",
                "05_screen_cdr3_cross_patient_similarity.py",
                "py",
                "Screens probed CDR3 sequences for cross-patient similarity (public/quasi-public TCR motifs) that could cause probe cross-reactivity.",
                "reports/tcr/cdr3_similarity_screen.tsv",
            ),
            (
                "08.06",
                "06_resolve_multiclonal_and_ambiguous_cells.py",
                "py",
                "Classifies singlet, probable multiplet, low-confidence and unassigned TCR calls.",
                "data/derived/tcr_resolved_calls.parquet",
            ),
            (
                "08.07",
                "07_build_clone_metadata_table.py",
                "py",
                "Summarises clone size, patient, section support, phenotype composition and replicate recurrence.",
                "data/derived/clone_metadata.parquet",
            ),
            (
                "08.08",
                "08_generate_tcr_release_report.py",
                "py",
                "Freezes high-confidence clone definitions for primary analysis and documents excluded/ambiguous calls.",
                "data/releases/v1_tcr_calls/",
            ),
        ],
    ),
    (
        9,
        "09_spatial_graph_construction_and_calibration",
        "Spatial graph construction, sensitivity framework and null-model calibration",
        [
            (
                "09.00",
                "00_generate_candidate_spatial_graphs.py",
                "py",
                "Builds radius, k-nearest-neighbour, Delaunay and boundary-contact graphs for every section.",
                "data/graphs/candidate/",
            ),
            (
                "09.01",
                "01_prune_graphs_for_tissue_gaps.py",
                "py",
                "Removes Delaunay/radius edges that bridge tissue folds, holes or QC-excluded regions.",
                "data/graphs/pruned/",
            ),
            (
                "09.02",
                "02_calibrate_graph_parameters.py",
                "py",
                "Selects biologically interpretable scales using cell size, interaction ranges, graph connectivity and stability.",
                "config/graph_parameters.yaml",
            ),
            (
                "09.03",
                "03_construct_primary_cell_graph.py",
                "py",
                "Creates the prespecified primary graph with node/edge metadata and patient-separated components.",
                "data/graphs/primary_graphs/",
            ),
            (
                "09.04",
                "04_construct_tumour_tcell_bipartite_graph.py",
                "py",
                "Represents direct and near-direct malignant-cell/T-cell relationships for contact-focused analyses.",
                "data/graphs/tumour_tcell/",
            ),
            (
                "09.05",
                "05_construct_clone_induced_subgraphs.py",
                "py",
                "Extracts all cells belonging to each clone plus local microenvironment shells.",
                "data/graphs/clones/",
            ),
            (
                "09.06",
                "06_run_graph_sensitivity_grid.py",
                "py",
                "Repeats core metrics over plausible radii/k values to identify robust conclusions.",
                "reports/graphs/sensitivity_grid.parquet",
            ),
            (
                "09.07",
                "07_generate_synthetic_ground_truth_patterns.py",
                "py",
                "Simulates spatial point patterns with known, controllable effect sizes and known clone-spatial relationships.",
                "data/graphs/synthetic/",
            ),
            (
                "09.08",
                "08_run_calibration_suite_on_synthetic_data.py",
                "py",
                "Confirms each planned null model achieves nominal Type I error and adequate power at n=10-patient scale before it is applied to observed data.",
                "reports/graphs/null_model_calibration.pdf",
            ),
        ],
    ),
    (
        10,
        "10_niche_and_ecosystem_discovery",
        "Discovery of spatial ecosystems and tissue domains",
        [
            (
                "10.00",
                "00_compute_cell_type_neighbourhood_enrichment.py",
                "py",
                "Tests pairwise cell-type adjacency using within-section constrained permutations (calibrated in Spatial Graph Construction and Calibration).",
                "data/derived/neighbourhood_enrichment.parquet",
            ),
            (
                "10.01",
                "01_compute_local_neighbourhood_compositions.py",
                "py",
                "Creates cell-centred neighbourhood composition vectors at several spatial scales.",
                "data/derived/local_compositions.parquet",
            ),
            (
                "10.02",
                "02_discover_neighbourhood_archetypes.R",
                "R",
                "Identifies recurrent neighbourhoods with consensus clustering or topic modelling and quantifies stability.",
                "data/derived/neighbourhood_archetypes.parquet",
            ),
            (
                "10.03",
                "03_segment_tissue_domains.py",
                "py",
                "Forms contiguous domains from neighbourhood states while preserving boundaries and rare niches.",
                "data/derived/tissue_domains.parquet",
            ),
            (
                "10.04",
                "04_annotate_ecosystems_with_blinded_rules.py",
                "py",
                "Assigns descriptive ecosystem labels only after unsupervised discovery, using a documented rubric.",
                "metadata/ecosystem_annotation.tsv",
            ),
            (
                "10.05",
                "05_quantify_ecosystem_abundance_and_topology.py",
                "py",
                "Measures abundance, fragmentation, interface length, mixing, compactness and relation to tumour borders.",
                "data/derived/ecosystem_metrics.parquet",
            ),
            (
                "10.06",
                "06_test_patient_recurrence.R",
                "R",
                "Fits patient-level models with technical replicates nested within patient; HPV is not tested here.",
                "reports/niches/recurrence_models.pdf",
            ),
            (
                "10.07",
                "07_leave_one_patient_out_niche_stability.py",
                "py",
                "Tests whether archetypes and labels remain identifiable when each patient is withheld.",
                "reports/niches/LOPO_stability.pdf",
            ),
        ],
    ),
    (
        11,
        "11_clone_spatial_descriptors",
        "Clone-level ecological descriptors and structure test (provisional)",
        [
            (
                "11.00",
                "00_compute_clone_spatial_descriptors_rarefied.py",
                "py",
                "Calculates dispersion, clustering, tumour distance, border enrichment and domain occupancy using rarefaction/size-normalisation.",
                "data/derived/clone_spatial_descriptors.parquet",
            ),
            (
                "11.01",
                "01_compute_clone_cell_state_composition.py",
                "py",
                "Summarises cytotoxic, exhausted, activated, proliferative and ambiguous states per clone with uncertainty.",
                "data/derived/clone_state_composition.parquet",
            ),
            (
                "11.02",
                "02_quantify_clone_tumour_engagement.py",
                "py",
                "Measures malignant-cell adjacency, contact normalised by opportunity, penetration and interface localisation; cross-checked against Quality Control spillover-risk and resegmentation outputs.",
                "data/derived/clone_tumour_engagement.parquet",
            ),
            (
                "11.03",
                "03_quantify_clone_apc_support.py",
                "py",
                "Measures proximity and neighbourhood enrichment with dendritic cells, macrophages and antigen-presentation programs.",
                "data/derived/clone_apc_support.parquet",
            ),
            (
                "11.04",
                "04_quantify_stromal_and_myeloid_barriers.py",
                "py",
                "Builds clone-specific fibroblast, vascular and suppressive-myeloid barrier metrics along paths to tumour, using calibrated null models.",
                "data/derived/clone_barrier_metrics.parquet",
            ),
            (
                "11.05",
                "05_test_discrete_vs_continuous_structure.R",
                "R",
                "Compares a discrete clustering model against a continuous archetypal/latent-trajectory model; mandatory before any taxonomy step.",
                "reports/clone_ecology/structure_test.pdf",
            ),
            (
                "11.06",
                "06_discover_provisional_structure.R",
                "R",
                "Fits whichever structure 11.05 supports (discrete fates or a continuous score), without imposing categories the data don't show.",
                "data/derived/clone_ecological_structure.parquet",
            ),
            (
                "11.07",
                "07_freeze_provisional_taxonomy_version.py",
                "py",
                "Freezes taxonomy/structure definitions as taxonomy_version=v1 and routes to External Checkpoint Validation for an external sanity check.",
                "governance/taxonomy_version_log.tsv",
            ),
        ],
    ),
    (
        12,
        "12_external_checkpoint_validation",
        "Early external sanity-check and taxonomy-freeze gate",
        [
            (
                "12.00",
                "00_project_provisional_signatures_to_bulk_reference.py",
                "py",
                "Projects the provisional clone-state/structure signatures onto an independent public HNSCC bulk or single-cell reference.",
                "reports/external_checkpoint/bulk_projection.pdf",
            ),
            (
                "12.01",
                "01_test_transcriptional_program_transfer.py",
                "py",
                "Tests whether the transcriptional programs distinguishing provisional structure categories are recoverable in the independent reference.",
                "reports/external_checkpoint/program_transfer.pdf",
            ),
            (
                "12.02",
                "02_quantify_directional_consistency.py",
                "py",
                "Quantifies whether the direction of any effect is consistent with the independent reference, not just whether p<0.05.",
                "results/external_checkpoint/directional_consistency.tsv",
            ),
            (
                "12.03",
                "03_decide_freeze_or_revise.py",
                "py",
                "Applies a predeclared decision rule: freeze taxonomy_version v1 and proceed, or revise to v2 (capped at one revision).",
                "governance/freeze_decision.tsv",
            ),
        ],
    ),
    (
        13,
        "13_clone_ecology_confirmatory_models",
        "Confirmatory variance-partitioning and clone ecological modelling",
        [
            (
                "13.00",
                "00_load_frozen_taxonomy_version.py",
                "py",
                "Loads the exact frozen taxonomy_version from External Checkpoint Validation; refuses to run if any upstream input hash has changed.",
                "data/releases/v1_clone_structure/",
            ),
            (
                "13.01",
                "01_fit_variance_partition_models.R",
                "R",
                "Fits hierarchical/Bayesian models decomposing variance in clone spatial/phenotypic descriptors into identity, context and patient components.",
                "reports/clone_ecology/variance_partition_models.pdf",
            ),
            (
                "13.02",
                "02_fit_hierarchical_clone_models.R",
                "R",
                "Fits mixed or Bayesian models for any confirmed categorical structure, nested in patient and section.",
                "reports/clone_ecology/hierarchical_models.pdf",
            ),
            (
                "13.03",
                "03_test_clone_size_as_confounder.R",
                "R",
                "Determines which associations survive after matching or adjustment for clone size and detection power.",
                "reports/clone_ecology/clone_size_sensitivity.pdf",
            ),
            (
                "13.04",
                "04_run_leave_one_patient_out_stability.R",
                "R",
                "Confirms the variance partition and any categorical structure are not driven by a single patient.",
                "reports/clone_ecology/lopo_stability.pdf",
            ),
            (
                "13.05",
                "05_test_segmentation_robustness.py",
                "py",
                "Re-runs the confirmatory models on the Quality Control resegmented reference subset and reports whether conclusions hold.",
                "reports/clone_ecology/segmentation_robustness.pdf",
            ),
            (
                "13.06",
                "06_generate_clone_atlas.py",
                "py",
                "Produces standardised maps and summaries for every high-confidence clone, preventing example-selection bias.",
                "reports/clone_ecology/clone_atlas.html",
            ),
        ],
    ),
    (
        14,
        "14_spatial_interactions_and_barriers",
        "Barrier topology and candidate spatial interactions",
        [
            (
                "14.00",
                "00_define_sender_receiver_pairs.py",
                "py",
                "Predeclares biologically motivated tumour-T cell, fibroblast-T cell, myeloid-T cell and APC-T cell comparisons.",
                "config/sender_receiver_pairs.yaml",
            ),
            (
                "14.01",
                "01_filter_ligand_receptor_database_to_panel.py",
                "py",
                "Restricts interaction databases to assayed genes and labels incomplete pairs to prevent overinterpretation.",
                "references/panel_supported_interactions.tsv",
            ),
            (
                "14.02",
                "02_compute_spatially_constrained_scores.py",
                "py",
                "Scores candidate-interaction support only across graph-connected cells and compares against degree-preserving nulls.",
                "data/derived/spatial_interaction_scores.parquet",
            ),
            (
                "14.03",
                "03_model_barrier_topology_by_structure.R",
                "R",
                "Tests whether the confirmed clone structure associates with checkpoint, chemokine, TGF-beta, interferon and antigen-presentation programs at barrier interfaces.",
                "reports/interactions/barrier_topology_models.pdf",
            ),
            (
                "14.04",
                "04_analyse_barrier_pathways.py",
                "py",
                "Evaluates programs specifically at spatial interfaces separating excluded clones from tumour.",
                "data/derived/barrier_pathways.parquet",
            ),
            (
                "14.05",
                "05_prioritise_testable_interactions.py",
                "py",
                "Ranks candidate interactions by effect size, spatial specificity, cross-patient consistency, panel completeness and external support. Never uses the word 'mechanism'.",
                "results/interaction_priority_table.tsv",
            ),
            (
                "14.06",
                "06_benchmark_against_published_barrier_studies.R",
                "R",
                "Quantitatively compares barrier-topology effect sizes against Grout/Kim 2022 and Hwang et al. 2022.",
                "reports/interactions/literature_benchmark.pdf",
            ),
        ],
    ),
    (
        15,
        "15_hpv_stratified_analysis",
        "HPV-stratified analysis (consolidated, capped)",
        [
            (
                "15.00",
                "00_validate_hpv_metadata_and_probe_signal.py",
                "py",
                "Cross-checks clinical HPV labels with oncoprotein probes and flags discordant or uncertain samples.",
                "metadata/hpv_status_validated.tsv",
            ),
            (
                "15.01",
                "01_prespecify_primary_hpv_contrasts.py",
                "py",
                "Fixes, before any test is run, a maximum of one or two primary HPV contrasts for the entire thesis.",
                "governance/hpv_primary_contrasts.yaml",
            ),
            (
                "15.02",
                "02_run_prospective_power_simulation.R",
                "R",
                "Simulates the minimum detectable effect size for n=5 vs n=5 under the planned model.",
                "reports/hpv/power_simulation.pdf",
            ),
            (
                "15.03",
                "03_compare_cellular_composition_patient_level.R",
                "R",
                "Compares HPV groups using patient as the sampling unit, not cells.",
                "reports/hpv/composition_models.pdf",
            ),
            (
                "15.04",
                "04_compare_ecosystem_and_clone_structure.R",
                "R",
                "Tests the single pre-specified contrast(s) against ecosystem abundance/topology and confirmed clone structure.",
                "reports/hpv/structure_models.pdf",
            ),
            (
                "15.05",
                "05_run_small_sample_robustness_checks.R",
                "R",
                "Uses permutation, exact tests, bootstrap and leave-one-patient-out analyses to expose fragility.",
                "reports/hpv/robustness.pdf",
            ),
            (
                "15.06",
                "06_prepare_hpv_claim_strength_table.py",
                "py",
                "Grades each HPV conclusion as exploratory, supported or unsuitable for inference.",
                "results/hpv_claim_strength.tsv",
            ),
        ],
    ),
    (
        16,
        "16_external_validation_and_generalisation",
        "Full external validation, triangulation and generalisability",
        [
            (
                "16.00",
                "00_define_validation_claims.py",
                "py",
                "Maps each major claim to a specific independent validation dataset and success criterion, predeclared before analysis.",
                "governance/validation_plan.tsv",
            ),
            (
                "16.01",
                "01_acquire_independent_spatial_dataset.py",
                "py",
                "Seeks at least one independent spatial dataset with clonal or lineage information to test whether the framework generalises.",
                "data/external/spatial/",
            ),
            (
                "16.02",
                "02_acquire_hnscc_scrna_references.py",
                "py",
                "Downloads selected public HNSCC scRNA-seq datasets with licences, metadata and checksums.",
                "data/external/scrna/",
            ),
            (
                "16.03",
                "03_validate_cell_state_signatures.R",
                "R",
                "Tests whether targeted-panel programs correspond to full-transcriptome T-cell and stromal states.",
                "reports/validation/signature_validation.pdf",
            ),
            (
                "16.04",
                "04_validate_ecosystem_signatures_in_bulk.py",
                "py",
                "Projects ecosystem-derived signatures to TCGA/CPTAC-like cohorts and cautiously tests clinical/immune associations.",
                "reports/validation/bulk_projection.pdf",
            ),
            (
                "16.05",
                "05_validate_framework_on_independent_dataset.py",
                "py",
                "Applies the full software framework end-to-end to the independent spatial dataset - the strongest single test of Q1.",
                "reports/validation/framework_generalisation.pdf",
            ),
            (
                "16.06",
                "06_compare_with_source_paper_results.py",
                "py",
                "Reproduces selected published McCord et al. results as a benchmark, then identifies which findings in this thesis are new.",
                "reports/validation/source_paper_comparison.tsv",
            ),
            (
                "16.07",
                "07_generate_evidence_matrix.py",
                "py",
                "Links every manuscript claim to discovery, sensitivity, replicate and external-validation evidence.",
                "results/claim_evidence_matrix.tsv",
            ),
        ],
    ),
    (
        17,
        "17_statistical_closure_and_release",
        "Statistical closure, figures, manuscript and reproducible release",
        [
            (
                "17.00",
                "00_freeze_primary_results.py",
                "py",
                "Locks the prespecified primary results (Q1-Q3, plus the single HPV contrast set) before exploratory extensions.",
                "data/releases/final_primary/",
            ),
            (
                "17.01",
                "01_control_multiplicity_and_report_effects.R",
                "R",
                "Applies a single, project-wide FDR/gatekeeping strategy across the prespecified primary family.",
                "results/statistical_summary.tsv",
            ),
            (
                "17.02",
                "02_generate_all_main_figures.py",
                "py",
                "Generates manuscript figures from frozen data with no manual editing of quantitative content.",
                "figures/main/",
            ),
            (
                "17.03",
                "03_generate_all_supplementary_figures.py",
                "py",
                "Creates complete diagnostics, sensitivity analyses, segmentation-robustness comparisons and non-selected clone maps.",
                "figures/supplementary/",
            ),
            (
                "17.04",
                "04_generate_results_tables.py",
                "py",
                "Creates publication-ready sample, QC, model and validation tables.",
                "tables/",
            ),
            (
                "17.05",
                "05_render_reproducible_manuscript.qmd",
                "qmd",
                "Builds the manuscript and methods supplement with embedded software and data provenance.",
                "manuscript/manuscript.pdf; methods.pdf",
            ),
            (
                "17.06",
                "06_build_public_data_release.py",
                "py",
                "Exports non-sensitive derived data, schemas and lightweight examples with licences.",
                "release/data/",
            ),
            (
                "17.07",
                "07_build_documented_software_package.py",
                "py",
                "Packages reusable functions, API documentation, tutorials and tests as an installable, versioned package.",
                "release/software/",
            ),
            (
                "17.08",
                "08_run_end_to_end_reproduction_test.sh",
                "sh",
                "Executes the full pipeline inside the container image, on a clean environment, first on a test subset and then the full dataset.",
                "reports/reproduction_test/",
            ),
            (
                "17.09",
                "09_run_null_model_calibration_regression.py",
                "py",
                "Re-runs the Spatial Graph Construction and Calibration simulation-calibration suite against the frozen release code, guarding against calibration drift.",
                "reports/release/calibration_regression.pdf",
            ),
            (
                "17.10",
                "10_archive_code_and_create_doi.sh",
                "sh",
                "Creates a tagged release, citation file, software bill of materials, container image digest and archive-ready bundle.",
                "release/archive_bundle/",
            ),
        ],
    ),
]

# Checkpoint grouping, following the project's "Core Snakemake execution
# sequence" specification -- some checkpoints deliberately span two
# phases (annotation_release = preprocess + annotation; tumour_and_tcr_release
# = tumour + tcr), matching that execution sequence exactly.
CHECKPOINTS = [
    ("project_ready", [1]),
    ("raw_data_frozen", [2]),
    ("canonical_objects", [3]),
    ("qc_and_integrity_release", [4]),
    ("annotation_release", [5, 6]),
    ("tumour_and_tcr_release", [7, 8]),
    ("graph_and_calibration_release", [9]),
    ("niche_release", [10]),
    ("clone_descriptor_release", [11]),
    ("external_checkpoint", [12]),
    ("clone_model_release", [13]),
    ("interaction_release", [14]),
    ("hpv_release", [15]),
    ("full_validation_release", [16]),
    ("publication_release", [17]),
]

SHELL_BY_LANG = {
    "py": "python3 {script}",
    "R": "Rscript {script}",
    "sh": "bash {script}",
    "qmd": "quarto render {script}",
}

TOP_LEVEL_DIRS = [
    "src",
    "workflow/rules",
    "config",
    "manifests",
    "governance",
    "environment/conda",
    "containers",
    "references",
    "tests/unit",
    "tests/simulation",
    "tests/fixtures/tiny_xenium",
    "notebooks/exploratory_only",
    "data/raw",
    "data/staged",
    "data/objects",
    "data/derived",
    "data/graphs",
    "data/releases",
    "data/external",
    "results/tables",
    "results/figures",
    "results/logs",
    "reports",
    "figures/main",
    "figures/supplementary",
    "tables",
    "manuscript/ch1_framework",
    "manuscript/ch2_variance_partition",
    "manuscript/ch3_barrier_topology",
    "manuscript/ch4_validation",
    "docs",
    "release",
]

PY_STUB_TEMPLATE = '''#!/usr/bin/env python3
"""
Phase {phase_idx:02d}.{order_suffix} -- {filename}

{purpose}

Primary output: {output}

Status: scaffold only -- no analytical logic implemented yet. Do not
implement analytical logic here without first checking
governance/analysis_registry.tsv for this script's registered
hypothesis, unit of analysis, and multiplicity family.
"""

from __future__ import annotations

import argparse
import sys

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.paths import find_project_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--project-root", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {{exc}}", file=sys.stderr)
        return 1

    raise NotImplementedError(
        "Phase {order_code} ({filename}) is not yet implemented."
    )


if __name__ == "__main__":
    sys.exit(main())
'''

R_STUB_TEMPLATE = """#!/usr/bin/env Rscript
# Phase {phase_idx:02d}.{order_suffix} -- {filename}
#
# {purpose}
#
# Primary output: {output}
#
# Status: scaffold only -- no analytical logic implemented yet.

stop("Phase {order_code} ({filename}) is not yet implemented.")
"""

SH_STUB_TEMPLATE = """#!/usr/bin/env bash
# Phase {phase_idx:02d}.{order_suffix} -- {filename}
#
# {purpose}
#
# Primary output: {output}
#
# Status: scaffold only -- no logic implemented yet.
set -Eeuo pipefail
echo "[ERROR] Phase {order_code} ({filename}) is not yet implemented." >&2
exit 1
"""

QMD_STUB_TEMPLATE = """---
title: "Phase {phase_idx:02d}.{order_suffix} -- {filename}"
format: pdf
---

<!--
{purpose}

Primary output: {output}

Status: scaffold only.
-->

::: {{.callout-warning}}
This manuscript-rendering document is not yet implemented.
:::
"""

TEMPLATES = {
    "py": PY_STUB_TEMPLATE,
    "R": R_STUB_TEMPLATE,
    "sh": SH_STUB_TEMPLATE,
    "qmd": QMD_STUB_TEMPLATE,
}


def write_file(path: Path, content: str, force: bool) -> str:
    if path.is_file() and not force:
        return "skip"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    if path.suffix in (".sh", ".py", ".R"):
        path.chmod(path.stat().st_mode | 0o111)
    return "wrote"


def _sentinel(folder: str, filename: str) -> str:
    stem = Path(filename).stem
    return f"results/logs/{folder}/.sentinels/{stem}.done"


def generate_workflow_rules(project_root: Path, force: bool) -> int:
    """Generate one .smk file per phase (sentinel-output rules, one per
    script) plus a _checkpoints.smk chaining the 16 named checkpoints in
    execution-sequence order, and a top-level Snakefile that includes all
    of it.

    Deliberately uses a sentinel-file pattern (touch a .done marker) rather
    than fabricating input/output wildcards for outputs the specification
    describes as directories or per-section paths (e.g. <section_id>.zarr):
    inventing that wiring now would be scientific/engineering work belonging
    to each phase's implementation, not to a scaffold. This still yields a
    valid, runnable Snakemake DAG end to end today.
    """
    written = 0
    rules_dir = project_root / "workflow" / "rules"
    rules_dir.mkdir(parents=True, exist_ok=True)

    phase_by_idx = {idx: (folder, purpose, scripts) for idx, folder, purpose, scripts in PHASES}

    for phase_idx, folder, purpose, scripts in PHASES:
        lines = [
            f'"""Phase {phase_idx:02d} ({folder}): {purpose}."""',
            "",
        ]
        for order_code, filename, lang, what, output in scripts:
            rule_name = f"phase{phase_idx:02d}_{Path(filename).stem}"
            script_rel = f"scripts/{folder}/{filename}"
            sentinel = _sentinel(folder, filename)

            shell_cmd = SHELL_BY_LANG[lang].format(script=script_rel)
            lines += [
                f"rule {rule_name}:",
                f'    """{order_code} -- {what} Primary output: {output}."""',
                f'    output: touch("{sentinel}")',
                f'    shell: "{shell_cmd}"',
                "",
                "",
            ]
        status = write_file(rules_dir / f"{folder}.smk", "\n".join(lines), force)
        written += status == "wrote"

    # Checkpoint chain: each checkpoint's rule requires every sentinel from
    # its constituent phase(s), and requires the previous checkpoint's own
    # output, so `snakemake <checkpoint>` runs everything in order.
    lines = [
        '"""Named checkpoints, matching the project\'s Core Snakemake execution sequence."""',
        "",
    ]
    previous_output = None
    for name, phase_indices in CHECKPOINTS:
        inputs = []
        if previous_output:
            inputs.append(f'"{previous_output}"')
        for pidx in phase_indices:
            folder, _purpose, scripts = (
                phase_by_idx[pidx][0],
                phase_by_idx[pidx][1],
                phase_by_idx[pidx][2],
            )
            for _order, filename, _lang, _what, _output in scripts:
                inputs.append(f'"{_sentinel(folder, filename)}"')
        checkpoint_output = f"results/logs/.checkpoints/{name}.done"
        input_block = ",\n        ".join(inputs)
        lines += [
            f"rule {name}:",
            "    input:",
            f"        {input_block}",
            f'    output: touch("{checkpoint_output}")',
            "",
            "",
        ]
        previous_output = checkpoint_output
    status = write_file(rules_dir / "_checkpoints.smk", "\n".join(lines), force)
    written += status == "wrote"

    return written


def scaffold(project_root: Path, force: bool, only: str = "all") -> None:
    """`only` scopes what --force is allowed to overwrite:
    - "all" (default): directories, script stubs, and workflow rules.
    - "rules": workflow/rules/*.smk only. Use this to regenerate rules
      (e.g. after editing CHECKPOINTS or fixing a rule-generation bug)
      without --force also clobbering script stubs that have since been
      implemented -- see git history around 2026-07-10 for why this flag
      exists: an unscoped --force wiped five Literature and Novelty Audit
      script implementations during that phase's milestone work.
    - "scripts": script stubs only, never touches workflow/rules/.
    """
    counts = {"wrote": 0, "skip": 0}

    if only not in ("all", "rules", "scripts"):
        raise ValueError(f"only must be 'all', 'rules', or 'scripts', got {only!r}")

    if only == "rules":
        n_rules = generate_workflow_rules(project_root, force)
        print(f"[OK]   Rules-only regeneration complete: {n_rules} rule file(s) touched.")
        return

    for d in TOP_LEVEL_DIRS:
        (project_root / d).mkdir(parents=True, exist_ok=True)

    for phase_idx, folder, purpose, scripts in PHASES:
        phase_dir = project_root / "scripts" / folder
        phase_dir.mkdir(parents=True, exist_ok=True)
        (project_root / "results" / "tables" / folder).mkdir(parents=True, exist_ok=True)
        (project_root / "results" / "figures" / folder).mkdir(parents=True, exist_ok=True)
        (project_root / "results" / "logs" / folder).mkdir(parents=True, exist_ok=True)

        for order_code, filename, lang, what, output in scripts:
            order_suffix = order_code.split(".")[1]
            content = TEMPLATES[lang].format(
                phase_idx=phase_idx,
                order_suffix=order_suffix,
                order_code=order_code,
                filename=filename,
                purpose=what,
                output=output,
                phase_name=purpose,
            )
            status = write_file(phase_dir / filename, content, force)
            counts[status] += 1

    # .gitkeep for every otherwise-empty data/results directory
    for base in [
        "data/raw",
        "data/staged",
        "data/objects",
        "data/derived",
        "data/graphs",
        "data/releases",
        "data/external",
        "notebooks/exploratory_only",
        "release",
    ]:
        write_file(project_root / base / ".gitkeep", "", force)

    if only == "all":
        n_rules = generate_workflow_rules(project_root, force)
        counts["wrote"] += n_rules

    print(
        f"[OK]   Scaffold complete: {counts['wrote']} files written, {counts['skip']} already present (kept)."
    )


def verify_policy(project_root: Path) -> int:
    """Minimal policy pass, mirroring the adopted legacy scaffold-as-linter
    pattern. Extended checks (taxonomy-freeze gate, HPV-contrast-count
    enforcement) belong in xenium_tcr_ecology.infra.policy once governance
    registries are populated with content, not at scaffold time."""
    errors = 0
    forbidden_path_fragments = ["/home/", "/root/", "/workspace/", "/mnt/"]
    scripts_dir = project_root / "scripts"
    for path in scripts_dir.rglob("*"):
        if not path.is_file():
            continue
        try:
            text = path.read_text()
        except UnicodeDecodeError:
            continue
        for frag in forbidden_path_fragments:
            if frag in text:
                print(f"[FAIL] Forbidden machine-specific path fragment '{frag}' in {path}")
                errors += 1
        if (
            "mechanism" in text.lower()
            and "14_spatial_interactions_and_barriers" not in str(path)
            and "07_tumour_epithelium_characterisation" not in str(path)
        ):
            print(
                f"[FAIL] Forbidden term 'mechanism' outside Spatial Interactions and Barriers/CNV appendix in {path}"
            )
            errors += 1

    required = [
        "src",
        "workflow/rules",
        "config",
        "manifests",
        "governance",
        "environment/conda",
        "containers",
        "tests",
        "scripts",
    ]
    for r in required:
        if not (project_root / r).is_dir():
            print(f"[FAIL] Missing required top-level directory: {r}")
            errors += 1

    if errors == 0:
        print("[OK]   Policy verification passed.")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--force", action="store_true", help="Overwrite existing files.")
    parser.add_argument(
        "--only",
        choices=["all", "rules", "scripts"],
        default="all",
        help="Scope --force to just workflow/rules (safe to use once scripts are implemented) "
        "or just script stubs. Default 'all' touches both -- be careful combining with --force "
        "once any script has been implemented.",
    )
    args = parser.parse_args()

    project_root = (
        Path(args.project_root).resolve()
        if args.project_root
        else Path(__file__).resolve().parents[1]
    )
    print(f"[INFO] Scaffolding at: {project_root}")

    scaffold(project_root, args.force, only=args.only)
    if args.only == "rules":
        return 0
    errors = verify_policy(project_root)
    if errors:
        print(f"[ERROR] Scaffold completed with {errors} policy violation(s).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
