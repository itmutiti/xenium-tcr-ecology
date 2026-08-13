"""Cross-checks clinical HPV/p16 labels against Xenium HPV16
oncoprotein probe signal (`15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py`) -- the gate `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`'s
prespecified contrast definition must pass through first.

**Clinically-grounded scope:** p16 immunohistochemistry is a validated
HPV surrogate marker only for oropharyngeal SCC (established clinical
practice, per CAP guidelines) -- this project's `sample_manifest.tsv`
reflects exactly this: `p16_ihc_status` is non-missing only for
oropharynx-site patients (P01, P09, P12, P13, P20, P28), and `NaN` for
the non-oropharyngeal sites (P15 Oral Cavity, P17 Oral Cavity, P19
Larynx, P23 Oral Cavity) where p16 IHC is not clinically indicated --
expected missingness, not a data-quality defect.

**Distinction that matters downstream:** the boolean `hpv_p16_positive`
field defaults to `False` for every one of these untested,
non-oropharyngeal patients -- a coding convention, not a confirmed-
negative determination. `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`'s contrast definition must not treat an
untested "presumed negative" with the same confidence as `P20`'s
clinically confirmed p16-negative result.

**Independent molecular check, where the Xenium panel allows it:** this
project's 623-gene panel includes 8 HPV16 probes (E1/E2/E4/E5/E6/E7/L1/
L2), but -- checked against `metadata/feature_annotation.tsv`'s
per-gene `sections_present` list, not assumed uniform -- these probes
are only present on 12/18 sections (missing from P15, P19 x2, P20 x2,
P23). `HPV16_E6`/`HPV16_E7` (the oncogenes whose continued expression
is mechanistically required for the malignant phenotype in HPV-driven
cancer, unlike the structural late genes L1/L2 which are frequently
lost on HPV genome integration) are used as the molecular check here,
not the full 8-probe panel.

**Data-driven threshold, not an arbitrary round number:** checking the
per-patient fraction of cells with detectable (`>0` raw count) E6/E7
signal reveals a clean natural gap -- {P09: 0.977, P12: 0.700, P13:
0.182, P17: 0.195} vs. {P28: 0.007, P10: 0.007, P01: 0.002} (true zeros
for P15/P19/P20/P23, which lack probe coverage) -- a ~26-fold separation
between the lowest "signal" patient (P13, 0.182) and the highest
"near-zero" patient (P28, 0.007). `PROBE_POSITIVE_FRACTION_THRESHOLD =
0.05` sits in this gap, not chosen blindly.

**Finding this validation surfaces:** two clinically p16-positive
patients (P01, P28) show near-zero HPV16 E6/E7 probe signal -- a
discordance flagged here, not silently inherited by `15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`'s contrast
definition. One clinically-untested patient (P17, Oral Cavity) shows a
moderate, unambiguous HPV16 E6/E7 signal (0.195) despite its
presumed-negative default label -- an unexpected positive molecular
signal in a patient p16 IHC was never clinically indicated for.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

HPV_ONCOGENE_PROBES = ["HPV16_E6", "HPV16_E7"]
PROBE_POSITIVE_FRACTION_THRESHOLD = 0.05


def compute_fraction_probe_positive(counts: np.ndarray) -> float:
    """Fraction of cells (rows) with a detectable (`>0`) count in any
    of the given gene columns."""
    return float((counts > 0).any(axis=1).mean())


def classify_validated_hpv_status(
    clinical_status: str,
    has_probe_coverage: bool,
    probe_positive_fraction: float | None,
    threshold: float = PROBE_POSITIVE_FRACTION_THRESHOLD,
) -> str:
    """Exhaustive 9-way classification combining the clinical p16 label
    (`"Positive"`/`"Negative"`/`"Not Tested"`), whether this patient's
    sections have any HPV16 probe coverage at all, and (when covered)
    whether the E6/E7 probe-positive fraction clears `threshold`."""
    is_probe_positive = (
        has_probe_coverage
        and probe_positive_fraction is not None
        and probe_positive_fraction >= threshold
    )

    if clinical_status == "Positive":
        if not has_probe_coverage:
            return "clinical_positive_no_molecular_verification"
        return (
            "confirmed_positive"
            if is_probe_positive
            else "discordant_clinical_positive_probe_negative"
        )
    if clinical_status == "Negative":
        if not has_probe_coverage:
            return "confirmed_negative_no_molecular_verification"
        return (
            "discordant_clinical_negative_probe_positive"
            if is_probe_positive
            else "confirmed_negative"
        )
    # clinical_status == "Not Tested"
    if not has_probe_coverage:
        return "presumed_negative_unverifiable"
    return (
        "probe_positive_clinically_untested"
        if is_probe_positive
        else "probe_negative_clinically_untested"
    )


def build_hpv_status_validation(project_root: Path) -> dict:
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    output_path = project_root / "metadata" / "hpv_status_validated.tsv"

    for path, phase in [
        (sample_manifest_path, None),
        (feature_annotation_path, None),
        (
            matrix_path,
            "`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found." + (f" Run {phase} first." if phase else ""))

    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
    primary = sample_manifest[sample_manifest["included_in_primary_hnscc_cohort"]].copy()
    primary["p16_ihc_status"] = primary["p16_ihc_status"].fillna("Not Tested")

    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")
    missing_probes = [
        g for g in HPV_ONCOGENE_PROBES if g not in set(feature_annotation["feature_name"])
    ]
    if missing_probes:
        raise PipelineError(
            f"HPV oncogene probe(s) {missing_probes} not found in '{feature_annotation_path}'."
        )
    probe_sections: set[str] = set()
    for gene in HPV_ONCOGENE_PROBES:
        row = feature_annotation.loc[feature_annotation["feature_name"] == gene].iloc[0]
        probe_sections |= set(str(row["sections_present"]).split(";"))

    adata = ad.read_h5ad(matrix_path)
    covered_cells = adata.obs["section_id"].isin(probe_sections)
    covered_subset = adata[covered_cells, HPV_ONCOGENE_PROBES]
    counts = covered_subset.layers["counts"]
    counts = counts.toarray() if hasattr(counts, "toarray") else counts
    cell_patient_ids = covered_subset.obs["patient_id"].to_numpy()

    per_patient_fraction: dict[str, float] = {}
    for patient_id in np.unique(cell_patient_ids):
        mask = cell_patient_ids == patient_id
        per_patient_fraction[patient_id] = compute_fraction_probe_positive(counts[mask])

    patient_probe_coverage = primary.groupby("patient_id")["section_id"].apply(
        lambda sections: bool(set(sections) & probe_sections)
    )
    patient_clinical_status = primary.groupby("patient_id")["p16_ihc_status"].first()

    rows = []
    for patient_id in patient_clinical_status.index:
        has_coverage = bool(patient_probe_coverage[patient_id])
        fraction = per_patient_fraction.get(patient_id) if has_coverage else None
        validated_status = classify_validated_hpv_status(
            patient_clinical_status[patient_id], has_coverage, fraction
        )
        rows.append(
            {
                "patient_id": patient_id,
                "clinical_p16_status": patient_clinical_status[patient_id],
                "has_hpv_probe_coverage": has_coverage,
                "hpv_e6_e7_probe_positive_fraction": fraction,
                "validated_hpv_status": validated_status,
            }
        )
    result = pd.DataFrame(rows).sort_values("patient_id").reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "n_patients": len(result),
        "n_confirmed_positive": int((result["validated_hpv_status"] == "confirmed_positive").sum()),
        "n_confirmed_negative_or_no_verification": int(
            result["validated_hpv_status"]
            .isin(["confirmed_negative", "confirmed_negative_no_molecular_verification"])
            .sum()
        ),
        "n_discordant": int(result["validated_hpv_status"].str.startswith("discordant").sum()),
        "n_presumed_negative_unverifiable": int(
            (result["validated_hpv_status"] == "presumed_negative_unverifiable").sum()
        ),
        "n_probe_positive_clinically_untested": int(
            (result["validated_hpv_status"] == "probe_positive_clinically_untested").sum()
        ),
        "output_path": str(output_path),
    }
