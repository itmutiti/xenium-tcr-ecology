"""Top-level Snakefile: includes one rule file per computational workflow
stage plus the named checkpoint chain.

This DAG covers raw data through final biological results only. The
`publication_release` target below (public data export and permanent
DOI/archive creation) sits outside it and is never part of the default
target: these actions are external-facing and, for the DOI, irreversible.

Run:
    snakemake --cores 8                        # full computational analysis (default target)
    snakemake --cores 8 <checkpoint_name>       # a specific named checkpoint
    snakemake --cores 8 publication_release     # public data export + code/DOI archive

See workflow/rules/*.smk for the stage-level rule definitions and
tools/scaffold_repository.py for how they are generated/regenerated.
"""

configfile: "config/config.yaml"

# Optional, additive: a single global container declaration so every
# rule *can* run inside containers/Apptainer.def's image via Snakemake's
# own `--use-apptainer` integration (profiles/apptainer/config.v8+.yaml),
# without editing any of the 152 individual rules below. Silently
# ignored -- confirmed directly from Snakemake's own documented
# behaviour -- unless `--use-apptainer`/`--software-deployment-method
# apptainer` is passed explicitly; native execution (the default,
# canonical route) is completely unaffected by this line's presence,
# including on hosts where `xenium-tcr-ecology.sif` does not exist.
# tools/run_with_apptainer.sh (wrapping the whole `snakemake` invocation
# in one `apptainer exec`) is the primary, recommended Apptainer route;
# this per-rule integration is a second, secondary option -- see
# docs/execution_manual/EXECUTION_MANUAL.md "Reproducibility".
container: "xenium-tcr-ecology.sif"

include: "workflow/rules/01_project_setup_and_governance.smk"
include: "workflow/rules/02_raw_data_ingestion.smk"
include: "workflow/rules/03_spatialdata_import.smk"
include: "workflow/rules/04_quality_control.smk"
include: "workflow/rules/05_preprocessing_and_normalisation.smk"
include: "workflow/rules/06_cell_type_annotation.smk"
include: "workflow/rules/07_tumour_epithelium_characterisation.smk"
include: "workflow/rules/08_tcr_clonal_analysis.smk"
include: "workflow/rules/09_spatial_graph_construction_and_calibration.smk"
include: "workflow/rules/10_niche_and_ecosystem_discovery.smk"
include: "workflow/rules/11_clone_spatial_descriptors.smk"
include: "workflow/rules/12_external_checkpoint_validation.smk"
include: "workflow/rules/13_clone_ecology_confirmatory_models.smk"
include: "workflow/rules/14_spatial_interactions_and_barriers.smk"
include: "workflow/rules/15_hpv_stratified_analysis.smk"
include: "workflow/rules/16_external_validation_and_generalisation.smk"
include: "workflow/rules/17_statistical_closure_and_release.smk"
include: "workflow/rules/_checkpoints.smk"

rule all:
    input:
        "results/logs/.checkpoints/computational_analysis_complete.done"
