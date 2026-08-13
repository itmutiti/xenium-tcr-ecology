"""Fixes the single primary HPV contrast for this project, before any
HPV hypothesis test runs (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`).

**Cap enforced here, not just described:** `config/
reproducibility_policy.yaml`'s `gates.hpv_single_contrast` and
`governance/analysis_registry.tsv`'s `hpv_primary_contrast` reserved
slot allow at most 1-2 confirmatory HPV contrasts overall. This module
defines exactly one -- a molecularly-validated positive-vs-negative
grouping -- and does not reserve or define a second, consistent with
this project's preference for minimal, well-justified scope.

**Grouping definition, built from `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`'s
validated status, not the raw clinical label:** `confirmed_positive`
and `probe_positive_clinically_untested` (`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`'s
finding that P17 shows unambiguous HPV16 E6/E7 molecular signal despite
never having been clinically tested -- p16 IHC is not clinically
indicated outside the oropharynx) are both grouped as "positive": the
direct Xenium oncogene-transcript readout is arguably a more specific
measure of transcriptionally active HPV than the p16 IHC surrogate.
`confirmed_negative`, `confirmed_negative_no_molecular_verification`,
`presumed_negative_unverifiable` and `probe_negative_clinically_
untested` are all "negative". The two discordant patients
(`discordant_clinical_positive_probe_negative`,
`discordant_clinical_negative_probe_positive`) and any patient with a
clinical label but no molecular verification opportunity in the
positive direction (`clinical_positive_no_molecular_verification`) are
"excluded" -- not confidently assignable to either group without further
work outside this milestone's scope, and left out of the primary
contrast rather than forced into a group.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError

POSITIVE_STATUSES = {"confirmed_positive", "probe_positive_clinically_untested"}
NEGATIVE_STATUSES = {
    "confirmed_negative",
    "confirmed_negative_no_molecular_verification",
    "presumed_negative_unverifiable",
    "probe_negative_clinically_untested",
}
CONTRAST_NAME = "hpv_validated_positive_vs_negative"


def assign_contrast_group(validated_hpv_status: str) -> str:
    """Pure, testable: maps every `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py` `validated_hpv_
    status` category to a contrast group -- exhaustive over the 9
    categories `classify_validated_hpv_status` can produce."""
    if validated_hpv_status in POSITIVE_STATUSES:
        return "positive"
    if validated_hpv_status in NEGATIVE_STATUSES:
        return "negative"
    return "excluded"


def build_hpv_contrast_config(hpv_status: pd.DataFrame, registered_date: str) -> dict:
    """Pure, testable: single contrast entry (the pre-existing
    `primary_contrasts` list schema in `governance/hpv_
    primary_contrasts.yaml`'s scaffold: `contrast_id`,
    `registered_date`, `model`, `minimum_detectable_effect`) from a
    (or synthetic-for-testing) `hpv_status_validated.tsv`-shaped
    DataFrame. `model` and `minimum_detectable_effect` are `None`
    placeholders here -- `model` is `15_hpv_stratified_analysis/03_compare_cellular_composition_patient_level.R`, `15_hpv_stratified_analysis/04_compare_ecosystem_and_clone_structure.R`'s job to fix per
    outcome domain, `minimum_detectable_effect` is Phase 15.02's job to
    compute."""
    groups = hpv_status.assign(
        contrast_group=hpv_status["validated_hpv_status"].apply(assign_contrast_group)
    )

    positive_ids = sorted(groups.loc[groups["contrast_group"] == "positive", "patient_id"])
    negative_ids = sorted(groups.loc[groups["contrast_group"] == "negative", "patient_id"])
    excluded_ids = sorted(groups.loc[groups["contrast_group"] == "excluded", "patient_id"])

    if len(positive_ids) == 0:
        raise PipelineError(
            "HPV-positive group is empty -- cannot define a contrast with zero patients on one side."
        )
    if len(negative_ids) == 0:
        raise PipelineError(
            "HPV-negative group is empty -- cannot define a contrast with zero patients on one side."
        )

    return {
        "contrast_id": CONTRAST_NAME,
        "registered_date": registered_date,
        "model": None,
        "minimum_detectable_effect": None,
        "positive_group": {"patient_ids": positive_ids, "criteria": sorted(POSITIVE_STATUSES)},
        "negative_group": {"patient_ids": negative_ids, "criteria": sorted(NEGATIVE_STATUSES)},
        "excluded_patients": {
            "patient_ids": excluded_ids,
            "reason": "discordant or unverifiable clinical/molecular HPV evidence",
        },
        "n_positive": len(positive_ids),
        "n_negative": len(negative_ids),
        "n_excluded": len(excluded_ids),
        "source": "metadata/hpv_status_validated.tsv (`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`)",
    }


MAX_PRIMARY_CONTRASTS = 2


def build_primary_hpv_contrasts(project_root: Path, registered_date: str = "2026-07-11") -> dict:
    hpv_status_path = project_root / "metadata" / "hpv_status_validated.tsv"
    output_path = project_root / "governance" / "hpv_primary_contrasts.yaml"

    if not hpv_status_path.exists():
        raise PipelineError(
            f"'{hpv_status_path}' not found. Run `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py` first."
        )

    existing_contrasts: list[dict] = []
    if output_path.exists():
        existing = yaml.safe_load(output_path.read_text()) or {}
        existing_contrasts = existing.get("primary_contrasts") or []

    already_registered = next(
        (c for c in existing_contrasts if c.get("contrast_id") == CONTRAST_NAME), None
    )
    if already_registered is None:
        if len(existing_contrasts) >= MAX_PRIMARY_CONTRASTS:
            raise PipelineError(
                f"'{output_path}' already has {len(existing_contrasts)} registered contrast(s) -- the cap "
                f"(config/reproducibility_policy.yaml gates.hpv_single_contrast, max {MAX_PRIMARY_CONTRASTS} for the "
                "entire thesis) does not permit registering another without an explicit, documented exception."
            )
        hpv_status = pd.read_csv(hpv_status_path, sep="\t")
        config = build_hpv_contrast_config(hpv_status, registered_date)
        existing_contrasts.append(config)
    else:
        config = already_registered

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "# Written by scripts/15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`), before any HPV\n"
        "# hypothesis test is run. A maximum of two entries are permitted here for this project;\n"
        "# every other HPV comparison anywhere in the pipeline is exploratory only and must not reference\n"
        "# hpv_status in a confirmatory model formula.\n"
        + yaml.dump({"primary_contrasts": existing_contrasts}, sort_keys=False)
    )

    return {
        "contrast_name": config["contrast_id"],
        "n_positive": config["n_positive"],
        "n_negative": config["n_negative"],
        "n_excluded": config["n_excluded"],
        "output_path": str(output_path),
    }
