"""Grades every HPV-related conclusion this phase could plausibly make
as `supported`, `exploratory`, or `unsuitable_for_inference`
(`15_hpv_stratified_analysis/06_prepare_hpv_claim_strength_table.py`) -- the final synthesis this deliberately small (n=4
vs n=4) design requires.

**Grading rule, not ad hoc per claim:** data-level findings
(`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`'s per-patient clinical-vs-molecular
discordances) and methodological/simulation facts (`15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`'s
minimum-detectable-effect-size result) are direct measurements or
verified simulation outputs, not underpowered group inferences -- graded
`supported`. Every group-comparison claim (all 25 categories tested
across `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`) is never graded `supported` -- none reached
BH-corrected significance, and `15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`'s power simulation
already establishes this n=4-vs-4 design cannot support a strong
confirmatory claim regardless of a nominally low p-value. Of those 25,
only the 3 `15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R` chose to stress-test can be graded
`exploratory` (directionally stable under leave-one-out and a bootstrap
CI excluding zero); the other 22 -- never selected as a preliminary
candidate worth stress-testing -- are `unsuitable_for_inference`
directly.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError

CONTRAST_ID = "hpv_validated_positive_vs_negative"


def grade_group_comparison(stress_test_info: dict | None) -> str:
    """Pure, testable: grading rule for a single group-comparison
    category. Never returns `"supported"` -- see module docstring for
    why. `exploratory` requires both a bootstrap CI excluding zero and
    leave-one-out direction stability (0 flips); anything else,
    including a category never stress-tested at all, is
    `unsuitable_for_inference`."""
    if stress_test_info is None:
        return "unsuitable_for_inference"
    bootstrap_excludes_zero = not (
        stress_test_info["bootstrap_ci_low"] <= 0 <= stress_test_info["bootstrap_ci_high"]
    )
    lopo_stable = stress_test_info["n_lopo_direction_flips"] == 0
    if bootstrap_excludes_zero and lopo_stable:
        return "exploratory"
    return "unsuitable_for_inference"


def build_hpv_claim_strength_table(project_root: Path) -> dict:
    hpv_status_path = project_root / "metadata" / "hpv_status_validated.tsv"
    contrasts_path = project_root / "governance" / "hpv_primary_contrasts.yaml"
    composition_path = (
        project_root / "data" / "derived" / "hpv_composition_comparison_results.parquet"
    )
    structure_path = project_root / "data" / "derived" / "hpv_structure_comparison_results.parquet"
    robustness_path = project_root / "data" / "derived" / "hpv_robustness_summary.parquet"
    output_path = project_root / "results" / "hpv_claim_strength.tsv"

    for path, phase in [
        (
            hpv_status_path,
            "`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`",
        ),
        (
            contrasts_path,
            "`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`, `15_hpv_stratified_analysis/02_run_prospective_power_simulation.R`",
        ),
        (
            composition_path,
            "`15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`",
        ),
        (structure_path, "`15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`"),
        (robustness_path, "`15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R`"),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    hpv_status = pd.read_csv(hpv_status_path, sep="\t")
    contrasts = yaml.safe_load(contrasts_path.read_text())
    contrast = next(
        (c for c in contrasts["primary_contrasts"] if c["contrast_id"] == CONTRAST_ID), None
    )
    if contrast is None:
        raise PipelineError(f"Contrast '{CONTRAST_ID}' not found in '{contrasts_path}'.")
    composition = pd.read_parquet(composition_path)
    structure = pd.read_parquet(structure_path)
    robustness = pd.read_parquet(robustness_path)

    claims = []

    for _, row in hpv_status[
        hpv_status["validated_hpv_status"].str.startswith("discordant")
    ].iterrows():
        claims.append(
            {
                "claim_id": f"hpv_discordance_{row['patient_id']}",
                "claim": f"Patient {row['patient_id']}'s clinical p16 status ({row['clinical_p16_status']}) is discordant with its Xenium HPV16 E6/E7 probe signal ({row['hpv_e6_e7_probe_positive_fraction']:.4f} of cells).",
                "inferential_status": "supported",
                "basis": "`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py` direct per-patient molecular measurement, not a group inference.",
            }
        )
    for _, row in hpv_status[
        hpv_status["validated_hpv_status"] == "probe_positive_clinically_untested"
    ].iterrows():
        claims.append(
            {
                "claim_id": f"hpv_probe_positive_untested_{row['patient_id']}",
                "claim": f"Patient {row['patient_id']} shows detectable HPV16 E6/E7 signal ({row['hpv_e6_e7_probe_positive_fraction']:.4f} of cells) despite no clinical p16 test.",
                "inferential_status": "supported",
                "basis": "`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py` direct per-patient molecular measurement, not a group inference.",
            }
        )

    mde = contrast.get("minimum_detectable_effect")
    claims.append(
        {
            "claim_id": "hpv_contrast_severely_underpowered",
            "claim": f"The primary HPV contrast (n={contrast['n_positive']} vs n={contrast['n_negative']}) requires an effect size of Cohen's d>={mde} for 80% power; an exact rank test at this n cannot reach p<0.05 without near-perfect group separation.",
            "inferential_status": "supported",
            "basis": "`15_hpv_stratified_analysis/02_run_prospective_power_simulation.R` verified Monte Carlo simulation and exact combinatorial calculation.",
        }
    )

    stress_test_lookup = {row["outcome"]: row.to_dict() for _, row in robustness.iterrows()}

    for _, row in composition.iterrows():
        stress_info = stress_test_lookup.get(row["lineage"])
        claims.append(
            {
                "claim_id": f"hpv_composition_{row['lineage']}",
                "claim": f"{row['lineage']} composition differs by HPV status (positive median={row['median_positive']:.4f}, negative median={row['median_negative']:.4f}, BH q={row['pvalue_bh']:.4f}).",
                "inferential_status": grade_group_comparison(stress_info),
                "basis": "`15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R` patient-level exact Wilcoxon test"
                + (
                    " + `15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R` stress test"
                    if stress_info
                    else ""
                ),
            }
        )

    for _, row in structure.iterrows():
        stress_key = (
            "clone_structure" if row["outcome_domain"] == "clone_structure" else row["category"]
        )
        stress_info = stress_test_lookup.get(stress_key)
        claims.append(
            {
                "claim_id": f"hpv_{row['outcome_domain']}_{row['category']}_{row['metric']}",
                "claim": f"{row['outcome_domain']} ({row['category']}, {row['metric']}) differs by HPV status (positive median={row['median_positive']:.4f}, negative median={row['median_negative']:.4f}, BH q={row['pvalue_bh']:.4f}).",
                "inferential_status": grade_group_comparison(stress_info),
                "basis": "`15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R` patient-level exact Wilcoxon test"
                + (
                    " + `15_hpv_stratified_analysis/05_run_small_sample_robustness_checks.R` stress test"
                    if stress_info
                    else ""
                ),
            }
        )

    result = pd.DataFrame(claims)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "n_claims": len(result),
        "n_supported": int((result["inferential_status"] == "supported").sum()),
        "n_exploratory": int((result["inferential_status"] == "exploratory").sum()),
        "n_unsuitable_for_inference": int(
            (result["inferential_status"] == "unsuitable_for_inference").sum()
        ),
        "output_path": str(output_path),
    }
