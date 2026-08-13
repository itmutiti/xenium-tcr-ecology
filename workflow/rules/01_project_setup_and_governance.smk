"""Project Setup and Governance (01_project_setup_and_governance): Research governance, computational environment and reproducibility."""

rule phase01_00_validate_project_scope:
    """01.00 -- Creates a machine-readable project charter recording the primary research questions; separates published findings from original hypotheses. Primary output: project_charter.yaml.
    Always succeeds -- see src/xenium_tcr_ecology/governance/charter.py."""
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/00_validate_project_scope.done")
    shell: "python3 scripts/01_project_setup_and_governance/00_validate_project_scope.py"


rule phase01_01_build_sample_manifest:
    """01.01 -- Parses GEO sample metadata into one canonical manifest with patient, run, replicate, HPV status, specimen type, file paths and checksums, and records the exact GEO record revision date used. Primary output: metadata/sample_manifest.tsv."""
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/01_build_sample_manifest.done")
    shell: "python3 scripts/01_project_setup_and_governance/01_build_sample_manifest.py"


rule phase01_02_generate_data_dictionary:
    """01.02 -- Defines every metadata field, unit, allowed value, missingness code and derivation rule. Primary output: metadata/data_dictionary.xlsx."""
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/02_generate_data_dictionary.done")
    shell: "python3 scripts/01_project_setup_and_governance/02_generate_data_dictionary.py"


rule phase01_03_lock_software_environments:
    """01.03 -- Builds pinned conda environments for Python and R; records OS and compiler information. Primary output: environment/conda/environment.resolved.yml; environment/conda/environment.lock; environment/R_sessionInfo.txt."""
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/03_lock_software_environments.done")
    shell: "bash scripts/01_project_setup_and_governance/03_lock_software_environments.sh"


rule phase01_04_build_container_images:
    """01.04 -- Builds a Docker/Apptainer image pinned to the locked environments, for the eventual frozen release. Primary output: containers/release.sif; containers/Dockerfile.
    Not required by `project_ready` (see that rule's own docstring in
    _checkpoints.smk): this rule's own script hard-fails without a
    reachable Docker daemon, which is not part of this project's
    documented base setup (README "Setup" covers conda + `.venv` only)
    and is confirmed absent on the Vast.ai instance types used for
    clean-room verification. Nothing in the computational DAG reads the
    built image. Run this by name on a host with Docker support
    when a container image is actually needed."""
    input: "results/logs/01_project_setup_and_governance/.sentinels/03_lock_software_environments.done"
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/04_build_container_images.done")
    shell: "bash scripts/01_project_setup_and_governance/04_build_container_images.sh"


rule phase01_05_initialise_reproducible_workflow:
    """01.05 -- Initialises Snakemake profiles, configuration templates, log directories, benchmark reporting and deterministic seeds. Primary output: Snakefile; config/config.yaml; profiles/."""
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/05_initialise_reproducible_workflow.done")
    shell: "python3 scripts/01_project_setup_and_governance/05_initialise_reproducible_workflow.py"


rule phase01_06_create_analysis_registry:
    """01.06 -- Registers each planned analysis with hypothesis, unit of analysis, exclusion criteria, primary endpoint, multiplicity family and taxonomy_version tracking. Primary output: governance/analysis_registry.tsv."""
    output: touch("results/logs/01_project_setup_and_governance/.sentinels/06_create_analysis_registry.done")
    shell: "python3 scripts/01_project_setup_and_governance/06_create_analysis_registry.py"
