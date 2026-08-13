"""Quality Control (04_quality_control): Multilevel quality control, principled exclusion and transcript-integrity robustness."""

rule phase04_00_compute_cell_level_qc_metrics:
    """04.00 -- Calculates transcript counts, detected genes, cell/nucleus area, control ratios and transcript density. Primary output: data/derived/cell_qc_metrics.parquet."""
    input: "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/04_quality_control/.sentinels/00_compute_cell_level_qc_metrics.done")
    shell: "python3 scripts/04_quality_control/00_compute_cell_level_qc_metrics.py"


rule phase04_01_compute_transcript_level_qc_metrics:
    """04.01 -- Profiles Q-values, unassigned transcripts, controls, nuclear overlap and spatial transcript-density artefacts. Primary output: data/derived/transcript_qc_metrics.parquet."""
    input: "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/04_quality_control/.sentinels/01_compute_transcript_level_qc_metrics.done")
    shell: "python3 scripts/04_quality_control/01_compute_transcript_level_qc_metrics.py"


rule phase04_02_detect_spatial_qc_artifacts:
    """04.02 -- Finds damaged regions, edge effects, holes, striping and local decoding failures using spatial statistics. Primary output: reports/qc/spatial_artifact_masks/."""
    input: "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/04_quality_control/.sentinels/02_detect_spatial_qc_artifacts.done")
    shell: "python3 scripts/04_quality_control/02_detect_spatial_qc_artifacts.py"


rule phase04_03_assess_segmentation_quality:
    """04.03 -- Measures cell geometry, nucleus containment, transcript assignment and suspicious multinucleation; samples cells for visual review. Primary output: reports/qc/segmentation_review.html."""
    input: "results/logs/.checkpoints/canonical_objects.done"
    output: touch("results/logs/04_quality_control/.sentinels/03_assess_segmentation_quality.done")
    shell: "python3 scripts/04_quality_control/03_assess_segmentation_quality.py"


rule phase04_04_estimate_transcript_spillover:
    """04.04 -- Flags cells near segmentation boundaries adjacent to a different predicted cell type; estimates a per-cell spillover-risk score using distance-to-boundary and neighbour-identity weighting. Primary output: data/derived/spillover_risk.parquet.
    Numbered under Quality Control (its topical phase) but runs after Cell Type Annotation:
    needs 06_cell_type_annotation/06_integrate_annotation_evidence.py's
    lineage calls.
    Deliberately excluded from qc_and_integrity_release's own input:
    (including it there would be cyclic -- this rule needs Cell Type Annotation, which
    itself needs qc_and_integrity_release -- and it neither gates
    the QC threshold/filter/release pipeline anyway; see that checkpoint
    rule's own docstring in workflow/rules/_checkpoints.smk for the full
    reasoning)."""
    input:
        "results/logs/.checkpoints/canonical_objects.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/04_quality_control/.sentinels/04_estimate_transcript_spillover.done")
    shell: "python3 scripts/04_quality_control/04_estimate_transcript_spillover.py"


rule phase04_05_resegment_reference_subset:
    """04.05 -- Independently resegments a representative subset of sections with an alternative transcript-reassignment method; used solely to test whether headline spatial results are robust to segmentation choice. Primary output: data/objects/resegmented_subset/.
    Deferred alongside 04.04 above as a thematic pairing, not because of its own data dependency on Cell Type Annotation -- checked directly
    against src/xenium_tcr_ecology/qc/resegmentation.py: it reads only
    canonical_objects/spatialdata and the primary analysis matrix, neither of
    which is a Cell Type Annotation output. Excluded from qc_and_integrity_release's
    input: for the same reason as 04.04 -- see that rule's docstring above
    and workflow/rules/_checkpoints.smk's qc_and_integrity_release docstring."""
    input:
        "results/logs/.checkpoints/canonical_objects.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done"
    output: touch("results/logs/04_quality_control/.sentinels/05_resegment_reference_subset.done")
    shell: "python3 scripts/04_quality_control/05_resegment_reference_subset.py"


rule phase04_06_define_qc_thresholds_hierarchically:
    """04.06 -- Uses robust, section-aware thresholds rather than global cut-offs; documents sensitivity alternatives. Primary output: config/qc_thresholds.yaml."""
    input:
        "results/logs/04_quality_control/.sentinels/00_compute_cell_level_qc_metrics.done",
        "results/logs/04_quality_control/.sentinels/02_detect_spatial_qc_artifacts.done"
    output: touch("results/logs/04_quality_control/.sentinels/06_define_qc_thresholds_hierarchically.done")
    shell: "Rscript scripts/04_quality_control/06_define_qc_thresholds_hierarchically.R"


rule phase04_07_apply_qc_filters_with_audit_trail:
    """04.07 -- Applies filters without deleting data and records a reason code for every excluded cell and transcript. Primary output: data/objects/qc_filtered.h5ad; data/derived/exclusion_log.tsv."""
    input:
        "results/logs/04_quality_control/.sentinels/00_compute_cell_level_qc_metrics.done",
        "results/logs/04_quality_control/.sentinels/02_detect_spatial_qc_artifacts.done",
        "results/logs/04_quality_control/.sentinels/06_define_qc_thresholds_hierarchically.done"
    output: touch("results/logs/04_quality_control/.sentinels/07_apply_qc_filters_with_audit_trail.done")
    shell: "python3 scripts/04_quality_control/07_apply_qc_filters_with_audit_trail.py"


rule phase04_08_assess_replicate_concordance:
    """04.08 -- Quantifies technical replicate agreement in counts, cell composition, gene expression and spatial statistics, correctly modelled as day nested within patient. Primary output: reports/qc/replicate_concordance.pdf."""
    input: "results/logs/04_quality_control/.sentinels/07_apply_qc_filters_with_audit_trail.done"
    output: touch("results/logs/04_quality_control/.sentinels/08_assess_replicate_concordance.done")
    shell: "Rscript scripts/04_quality_control/08_assess_replicate_concordance.R"


rule phase04_09_generate_qc_release_report:
    """04.09 -- Produces a section-by-section QC report and a go/no-go decision before biological analysis. Primary output: reports/qc/QC_release_report.html."""
    input:
        "results/logs/04_quality_control/.sentinels/00_compute_cell_level_qc_metrics.done",
        "results/logs/04_quality_control/.sentinels/01_compute_transcript_level_qc_metrics.done",
        "results/logs/04_quality_control/.sentinels/02_detect_spatial_qc_artifacts.done",
        "results/logs/04_quality_control/.sentinels/03_assess_segmentation_quality.done",
        "results/logs/04_quality_control/.sentinels/06_define_qc_thresholds_hierarchically.done",
        "results/logs/04_quality_control/.sentinels/07_apply_qc_filters_with_audit_trail.done",
        "results/logs/04_quality_control/.sentinels/08_assess_replicate_concordance.done"
    output: touch("results/logs/04_quality_control/.sentinels/09_generate_qc_release_report.done")
    shell: "python3 scripts/04_quality_control/09_generate_qc_release_report.py"
