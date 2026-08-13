"""TCR CDR3 probe registry construction (`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`).

Separates patient-specific CDR3 probes from conventional T-cell genes
using feature metadata (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s `classify_feature`) and the
probe naming convention (a date/batch prefix, the CDR3 amino acid
sequence, and the TCR chain -- confirmed against all 11 date/batch
prefixes in the panel),
then structures each probe into its parsed components plus the section(s)
/ patient(s) whose panel physically includes it
(`results/tables/03_spatialdata_import/gene_panel_membership.parquet`).

This module builds the registry; it does not audit it for leakage or
naming conflicts across patients (`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s explicit job) or assess
the ascertainment/selection process behind which clonotypes were probed
at all (`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`'s explicit job) -- both are separate, later milestones
per the project's phase breakdown, not pre-empted here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Same naming convention as CDR3_PROBE_PATTERN (feature_classification.py),
# with capture groups added for parsing -- kept as a distinct constant
# rather than modifying the shared classification pattern, since that
# pattern is also used (unparenthesised) elsewhere for a plain match test.
CDR3_PROBE_PATTERN_GROUPED = re.compile(r"^([0-9]{6}[A-Z]?)_([A-Z]+)_(TR[AB])$")


def parse_cdr3_probe_name(name: str) -> dict | None:
    """Pure, testable parse of one probe name into its structured
    components. Returns None for a name that does not match the CDR3
    probe convention (callers should already have restricted to
    `feature_class == 'tcr_cdr3_probe'` via `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`'s classification,
    so a non-match here indicates an inconsistency worth surfacing, not
    silently skipping)."""
    match = CDR3_PROBE_PATTERN_GROUPED.match(name)
    if match is None:
        return None
    date_batch_prefix, cdr3_amino_acid_sequence, tcr_chain = match.groups()
    return {
        "probe_name": name,
        "date_batch_prefix": date_batch_prefix,
        "cdr3_amino_acid_sequence": cdr3_amino_acid_sequence,
        "tcr_chain": tcr_chain,
    }


def build_tcr_probe_registry(project_root: Path) -> dict:
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    panel_membership_path = (
        project_root
        / "results"
        / "tables"
        / "03_spatialdata_import"
        / "gene_panel_membership.parquet"
    )
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    output_path = project_root / "metadata" / "tcr_probe_registry.tsv"

    for p in (feature_annotation_path, panel_membership_path, sample_manifest_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")
    cdr3_probes = feature_annotation.loc[
        feature_annotation["feature_class"] == "tcr_cdr3_probe", "feature_name"
    ]
    if len(cdr3_probes) == 0:
        raise PipelineError(
            "No 'tcr_cdr3_probe' features found in feature_annotation.tsv. Run `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py` first."
        )

    panel_membership = pd.read_parquet(panel_membership_path)
    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
    section_to_patient = sample_manifest.set_index("section_id")["patient_id"].to_dict()

    rows = []
    n_unparseable = 0
    for probe_name in sorted(cdr3_probes):
        parsed = parse_cdr3_probe_name(probe_name)
        if parsed is None:
            n_unparseable += 1
            continue

        if probe_name in panel_membership.index:
            sections_with_probe = panel_membership.columns[
                panel_membership.loc[probe_name]
            ].tolist()
        else:
            sections_with_probe = []
        patients_with_probe = sorted(
            {section_to_patient[s] for s in sections_with_probe if s in section_to_patient}
        )

        parsed["n_sections_with_probe"] = len(sections_with_probe)
        parsed["sections_with_probe"] = ";".join(sections_with_probe)
        parsed["n_patients_with_probe"] = len(patients_with_probe)
        parsed["patients_with_probe"] = ";".join(patients_with_probe)
        rows.append(parsed)

    if n_unparseable > 0:
        raise PipelineError(
            f"{n_unparseable} feature(s) classified as 'tcr_cdr3_probe' did not match the parsing pattern -- "
            "classification and parsing are inconsistent, investigate before proceeding."
        )

    registry = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    registry.to_csv(output_path, sep="\t", index=False)

    return {
        "n_probes": len(registry),
        "n_tra_probes": int((registry["tcr_chain"] == "TRA").sum()),
        "n_trb_probes": int((registry["tcr_chain"] == "TRB").sum()),
        "n_distinct_date_batch_prefixes": int(registry["date_batch_prefix"].nunique()),
        "n_probes_single_patient": int((registry["n_patients_with_probe"] == 1).sum()),
        "n_probes_multi_patient": int((registry["n_patients_with_probe"] > 1).sum()),
        "n_probes_zero_patients": int((registry["n_patients_with_probe"] == 0).sum()),
        "output_path": str(output_path),
    }
