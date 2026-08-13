"""Clone metadata table construction (`08_tcr_clonal_analysis/07_build_clone_metadata_table.py`).

Summarises clone size, patient, section support, phenotype composition
and replicate recurrence for every distinct clone identified in Phase
8.06's resolved calls.

**Clone identity, a judgment call made and documented directly:** a
"clone" is defined here as the exact set of detected probes
(`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s `detected_probes` string) among cells classified `singlet`
or `low_confidence` (`08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`) -- cells with a well-defined detected-
probe identity. `probable_multiplet` and `unassigned` cells are excluded
by construction (they have no single well-defined clone identity to
group by). This sidesteps needing to resolve whether a TRA+TRB pair
(`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s `likely_single_clone_tra_trb_pair`) truly represents one
paired receptor -- an unconfirmed limitation already documented in
`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`, `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py` -- rather than silently assuming a pairing scheme this
project has no data to verify.

**Internal-consistency check:** every cell contributing to
a given `detected_probes` value is, by `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s construction, restricted to
that probe's intended patient -- so every clone's `patient_id` should be
internally uniform by construction. Verified before building the table;
a violation would indicate a bug upstream,
not a normal data condition, and is treated as a hard error.

**Section support and replicate recurrence:** for each clone, the
distinct sections its cells appear in are recorded directly. For the 7
patients with a technical-replicate section pair (`metadata/sample_manifest.tsv`, already used for this exact
purpose in `04_quality_control/08_assess_replicate_concordance.R`'s replicate concordance analysis), `detected_in_both_replicates`
records whether the same clone (same detected-probe set) was independently
detected in both runs -- a direct reproducibility signal, not assumed to
hold.

**Phenotype composition:** joined from `06_cell_type_annotation/04_resolve_t_cell_substates.R`'s `t_cell_states.parquet`
(`t_cell_state` per cell) -- the fraction of each clone's cells in each
T-cell state.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError


def build_clone_metadata(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    t_cell_states_path = project_root / "data" / "derived" / "t_cell_states.parquet"
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    output_path = project_root / "data" / "derived" / "clone_metadata.parquet"

    for p in (resolved_calls_path, t_cell_states_path, sample_manifest_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    resolved = pd.read_parquet(resolved_calls_path)
    t_cell_states = pd.read_parquet(t_cell_states_path).set_index("cell_id")
    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")

    clonal = resolved[resolved["resolution"].isin(["singlet", "low_confidence"])].copy()
    if len(clonal) == 0:
        raise PipelineError(
            "No 'singlet' or 'low_confidence' cells found -- cannot build clone metadata."
        )

    clonal["clone_id"] = clonal["detected_probes"]
    clonal["t_cell_state"] = t_cell_states.reindex(clonal.index)["t_cell_state"]

    # Internal-consistency check (see module docstring) -- a clone
    # spanning more than one patient would indicate an upstream bug in
    # `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s patient-restriction logic, not a normal
    # condition to silently tolerate.
    patients_per_clone = clonal.groupby("clone_id")["patient_id"].nunique()
    if (patients_per_clone > 1).any():
        bad_clones = patients_per_clone[patients_per_clone > 1].index.tolist()
        raise PipelineError(
            f"{len(bad_clones)} clone(s) span more than one patient (e.g. {bad_clones[:3]}) -- "
            "this violates `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s patient-restriction invariant and indicates an upstream bug."
        )

    replicate_patients = set(
        sample_manifest.loc[sample_manifest["is_technical_replicate"], "patient_id"]
    )
    section_to_run = sample_manifest.set_index("section_id")["run_number"].to_dict()

    rows = []
    for clone_id, group in clonal.groupby("clone_id", observed=True):
        patient_id = group["patient_id"].iloc[0]
        sections = sorted(group["section_id"].unique())
        phenotype_counts = group["t_cell_state"].value_counts(normalize=True).round(4).to_dict()

        detected_in_both_replicates = None
        if patient_id in replicate_patients:
            runs_present = {section_to_run.get(s) for s in sections}
            detected_in_both_replicates = {1, 2}.issubset(runs_present)

        rows.append(
            {
                "clone_id": clone_id,
                "patient_id": patient_id,
                "n_cells": len(group),
                "n_singlet_cells": int((group["resolution"] == "singlet").sum()),
                "n_low_confidence_cells": int((group["resolution"] == "low_confidence").sum()),
                "n_sections": len(sections),
                "sections": ";".join(sections),
                "is_replicate_patient": patient_id in replicate_patients,
                "detected_in_both_replicates": detected_in_both_replicates,
                "phenotype_composition": phenotype_counts,
                "dominant_t_cell_state": (
                    max(phenotype_counts, key=phenotype_counts.get) if phenotype_counts else None
                ),
            }
        )

    clone_metadata = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    clone_metadata.to_parquet(output_path)

    replicate_clones = clone_metadata[clone_metadata["is_replicate_patient"]]
    return {
        "n_clones": len(clone_metadata),
        "n_cells_in_clones": int(clone_metadata["n_cells"].sum()),
        "median_clone_size": float(clone_metadata["n_cells"].median()),
        "max_clone_size": int(clone_metadata["n_cells"].max()),
        "n_clones_in_replicate_patients": len(replicate_clones),
        "n_clones_detected_in_both_replicates": int(
            replicate_clones["detected_in_both_replicates"].sum()
        ),
        "output_path": str(output_path),
    }
