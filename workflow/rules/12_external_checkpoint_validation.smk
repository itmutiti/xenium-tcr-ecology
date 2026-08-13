"""External Checkpoint Validation (12_external_checkpoint_validation): Early external sanity-check and taxonomy-freeze gate."""

rule phase12_00_project_provisional_signatures_to_bulk_reference:
    """12.00 -- Projects the provisional clone-state/structure signatures onto an independent public HNSCC bulk or single-cell reference. Primary output: reports/external_checkpoint/bulk_projection.pdf."""
    input: "results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done"
    output: touch("results/logs/12_external_checkpoint_validation/.sentinels/00_project_provisional_signatures_to_bulk_reference.done")
    shell: "python3 scripts/12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py"


rule phase12_01_test_transcriptional_program_transfer:
    """12.01 -- Tests whether the transcriptional programs distinguishing provisional structure categories are recoverable in the independent reference. Reads data/external/GSE103322/GSE103322_HNSCC_all_data.txt.gz directly (src/xenium_tcr_ecology/external_checkpoint/program_transfer.py's build_program_transfer_test), but does not acquire it itself -- GSE103322 acquisition is centralised in phase12_00 (the only caller of bulk_reference.ensure_gse103322_acquired), matching the same dependency phase12_02 and phase12_05 already declare below for the identical reason. Depending on phase12_00's sentinel here, rather than having this rule download the file independently, also avoids a second, unsafe hazard: infra/download.py's download_file() writes to a fixed '<dest>.part' temp path with no locking, so two rules downloading the same destination concurrently could corrupt each other's write, not just race on reading a partial file. Primary output: reports/external_checkpoint/program_transfer.pdf."""
    input:
        "results/logs/12_external_checkpoint_validation/.sentinels/00_project_provisional_signatures_to_bulk_reference.done",
        "results/logs/06_cell_type_annotation/.sentinels/04_resolve_t_cell_substates.done",
        "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done"
    output: touch("results/logs/12_external_checkpoint_validation/.sentinels/01_test_transcriptional_program_transfer.done")
    shell: "python3 scripts/12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py"


rule phase12_02_quantify_directional_consistency:
    """12.02 -- Quantifies whether the direction of any effect is consistent with the independent reference, not just whether p<0.05. Primary output: results/external_checkpoint/directional_consistency.tsv."""
    input:
        "results/logs/12_external_checkpoint_validation/.sentinels/00_project_provisional_signatures_to_bulk_reference.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/06_discover_provisional_structure.done"
    output: touch("results/logs/12_external_checkpoint_validation/.sentinels/02_quantify_directional_consistency.done")
    shell: "python3 scripts/12_external_checkpoint_validation/02_quantify_directional_consistency.py"


rule phase12_03_decide_freeze_or_revise:
    """12.03 -- Applies a predeclared decision rule: freeze taxonomy_version v1 and proceed, or revise to v2 (capped at one revision). Primary output: governance/freeze_decision.tsv."""
    input:
        "results/logs/12_external_checkpoint_validation/.sentinels/01_test_transcriptional_program_transfer.done",
        "results/logs/12_external_checkpoint_validation/.sentinels/02_quantify_directional_consistency.done",
        "results/logs/11_clone_spatial_descriptors/.sentinels/07_freeze_provisional_taxonomy_version.done"
    output: touch("results/logs/12_external_checkpoint_validation/.sentinels/03_decide_freeze_or_revise.done")
    shell: "python3 scripts/12_external_checkpoint_validation/03_decide_freeze_or_revise.py"


rule phase12_05_rescore_cycling_state_with_primary_method:
    """12.05 -- Post-hoc follow-up (added after further review, not part of the original 18-stage specification; see docs/analysis_amendments.md): re-scores GSE103322's T cells with the primary pipeline's scanpy.tl.score_genes method (rather than the original simplified z-score proxy) to test whether the 40.3%-vs-14.4% Cycling-state discrepancy is a scoring-method artefact. Primary output: reports/external_checkpoint/cycling_rescore_comparison.pdf."""
    input: "results/logs/12_external_checkpoint_validation/.sentinels/00_project_provisional_signatures_to_bulk_reference.done"
    output: touch("results/logs/12_external_checkpoint_validation/.sentinels/05_rescore_cycling_state_with_primary_method.done")
    shell: "python3 scripts/12_external_checkpoint_validation/05_rescore_cycling_state_with_primary_method.py"
