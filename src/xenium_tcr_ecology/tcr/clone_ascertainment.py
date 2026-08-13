"""Clone ascertainment documentation (`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`).

Records how each probed clonotype was selected relative to each patient's
full TCR repertoire, and publishes this as an explicit boundary condition
on every downstream generalisability claim -- per the blueprint's own
framing, this is fundamentally a documentation deliverable, not a fresh
statistical test.

**Real ascertainment criteria, taken directly from the source paper's own
stated methods** (McCord et al. 2026, Sci Immunol 11:eaec3133, not
invented here):
1. "Probe-based TCR detection prioritises abundant clonotypes and may
   miss rare functional clones" (Limitations/Discussion) -- the probed
   set is a deliberately ABUNDANCE-BIASED sample of each patient's
   repertoire, not a random or exhaustive one.
2. "VDJdb-matched microbial-reactive (EBV, influenza) TCRs are used as
   bystander-clone markers" (Figure 3) -- a minority of probes are
   deliberately-included bystander/viral-reactive controls, not all
   tumour-associated/expanded clones.
3. "Patient-specific TCR probes detect 61 validated clonotypes in >2,000
   T cells across 8 patients" (Figure 7) -- the source paper names and
   individually validates 61 clonotypes across 8 patients; this project's
   panel (216 probes across 11 patients, `08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`) is broader than that
   explicitly-named-and-validated set, so most of this project's own
   probes are not individually described in the source paper's text.
4. "Only modest correlation between scRNA-seq and Xenium-derived TCR
   frequencies" (Limitations) -- even for probes that are
   patient-matched, detected spatial abundance should not be read as a
   precise readout of true clonal frequency.

**An explicit data-availability limitation:** this project's only
companion dataset (GSE287301, `06_cell_type_annotation/03_map_external_scrna_reference.py`)
contains a gene-expression matrix only (`filtered_feature_bc_matrix`,
`patient_matrix.txt`, `aggregation.csv`) -- no VDJ/TCR-contig files exist
in that GEO deposit. There is therefore no way to independently
reconstruct any patient's full clonal repertoire or determine each probed
clone's true rank/selection statistic (e.g. "this was the 3rd most
expanded clone") from data available to this project. The ascertainment
record built here is consequently qualitative and cohort-level (the
criteria above, and `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s empirical detection-specificity
evidence), not a per-clone quantitative repertoire rank.

**Also explicitly out of scope here (left to the milestones that actually
own it):** whether a probed TRA sequence and a probed TRB sequence
represent one paired clone in the same cell is a cell-level detection
question (`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`, `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`, `08_tcr_clonal_analysis/07_build_clone_metadata_table.py`'s job, not resolved by probe-level
ascertainment documentation).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

ASCERTAINMENT_CRITERIA = (
    "Probe-based TCR detection prioritises abundant clonotypes from the companion scRNA/TCR-seq "
    "data and may miss rare functional clones; a minority of probes are deliberately-included "
    "VDJdb-matched bystander/viral-reactive (EBV, influenza) controls, not only tumour-associated "
    "expanded clones (McCord et al. 2026, Sci Immunol 11:eaec3133, Figure 3, Limitations). "
    "The probed set is therefore NOT a random or exhaustive sample of any patient's TCR repertoire."
)

DATA_AVAILABILITY_LIMITATION = (
    "This project's only companion dataset (GSE287301) contains a gene-expression matrix only, no "
    "VDJ/TCR-contig files -- confirmed against the GEO deposit's file listing. Each probed "
    "clone's true rank within its patient's full repertoire cannot be independently reconstructed or "
    "verified from data available to this project; this record is qualitative/cohort-level, not a "
    "per-clone quantitative repertoire rank."
)


def build_clone_ascertainment_record(project_root: Path) -> dict:
    registry_path = project_root / "metadata" / "tcr_probe_registry.tsv"
    audit_path = project_root / "reports" / "tcr" / "patient_probe_audit.tsv"
    output_path = project_root / "metadata" / "clone_ascertainment.tsv"

    for p in (registry_path, audit_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    registry = pd.read_csv(registry_path, sep="\t")
    audit = pd.read_csv(audit_path, sep="\t")

    record = registry.merge(
        audit[
            [
                "probe_name",
                "intended_patient",
                "intended_patient_identified",
                "top_patient_detection_rate",
            ]
        ],
        on="probe_name",
        how="left",
    )
    record["ascertainment_criteria"] = ASCERTAINMENT_CRITERIA
    record["data_availability_limitation"] = DATA_AVAILABILITY_LIMITATION
    record["individually_named_in_source_paper"] = (
        False  # see module docstring point 3; not independently verifiable per-probe
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    record.to_csv(output_path, sep="\t", index=False)

    return {
        "n_probes": len(record),
        "n_probes_with_empirically_identified_patient": int(
            record["intended_patient_identified"].sum()
        ),
        "n_probes_source_paper_named_validated_clonotypes": 61,
        "n_patients_source_paper_named": 8,
        "n_patients_this_project": int(
            registry["patients_with_probe"].str.split(";").explode().nunique()
        ),
        "output_path": str(output_path),
    }
