"""HPV-Stratified Analysis (15_hpv_stratified_analysis): HPV-stratified analysis (consolidated, capped)."""

rule phase15_00_validate_hpv_metadata_and_probe_signal:
    """15.00 -- Cross-checks clinical HPV labels with oncoprotein probes and flags discordant or uncertain samples. Primary output: metadata/hpv_status_validated.tsv."""
    input: "results/logs/05_preprocessing_and_normalisation/.sentinels/05_create_primary_analysis_matrix.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/00_validate_hpv_metadata_and_probe_signal.done")
    shell: "python3 scripts/15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py"


rule phase15_01_prespecify_primary_hpv_contrasts:
    """15.01 -- Fixes, before any test is run, a maximum of one or two primary HPV contrasts for the entire thesis. Primary output: governance/hpv_primary_contrasts.yaml."""
    input: "results/logs/15_hpv_stratified_analysis/.sentinels/00_validate_hpv_metadata_and_probe_signal.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done")
    shell: "python3 scripts/15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py"


rule phase15_02_run_prospective_power_simulation:
    """15.02 -- Simulates the minimum detectable effect size for n=5 vs n=5 under the planned model. Primary output: reports/hpv/power_simulation.pdf."""
    input: "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/02_run_prospective_power_simulation.done")
    shell: "Rscript scripts/15_hpv_stratified_analysis/02_run_prospective_power_simulation.R"


rule phase15_03_compare_cellular_composition_patient_level:
    """15.03 -- Compares HPV groups using patient as the sampling unit, not cells. Primary output: reports/hpv/composition_models.pdf."""
    input:
        "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/03_compare_cellular_composition_patient_level.done")
    shell: "Rscript scripts/15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R"


rule phase15_04_compare_ecosystem_and_clone_structure:
    """15.04 -- Tests the single pre-specified contrast(s) against ecosystem abundance/topology and confirmed clone structure. Primary output: reports/hpv/structure_models.pdf."""
    input:
        "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done",
        "results/logs/10_niche_and_ecosystem_discovery/.sentinels/05_quantify_ecosystem_abundance_and_topology.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/04_compare_ecosystem_and_clone_structure.done")
    shell: "Rscript scripts/15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R"


rule phase15_05_run_small_sample_robustness_checks:
    """15.05 -- Uses permutation, exact tests, bootstrap and leave-one-patient-out analyses to expose fragility. Primary output: reports/hpv/robustness.pdf."""
    input:
        "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done",
        "results/logs/06_cell_type_annotation/.sentinels/06_integrate_annotation_evidence.done",
        "results/logs/13_clone_ecology_confirmatory_models/.sentinels/00_load_frozen_taxonomy_version.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/05_run_small_sample_robustness_checks.done")
    shell: "Rscript scripts/15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R"


rule phase15_06_prepare_hpv_claim_strength_table:
    """15.06 -- Grades each HPV conclusion as exploratory, supported or unsuitable for inference. Primary output: results/hpv_claim_strength.tsv."""
    input:
        "results/logs/15_hpv_stratified_analysis/.sentinels/01_prespecify_primary_hpv_contrasts.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/03_compare_cellular_composition_patient_level.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/04_compare_ecosystem_and_clone_structure.done",
        "results/logs/15_hpv_stratified_analysis/.sentinels/05_run_small_sample_robustness_checks.done"
    output: touch("results/logs/15_hpv_stratified_analysis/.sentinels/06_prepare_hpv_claim_strength_table.done")
    shell: "python3 scripts/15_hpv_stratified_analysis/06_prepare_hpv_claim_strength_table.py"
