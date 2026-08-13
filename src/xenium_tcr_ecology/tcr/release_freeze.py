"""TCR calls freeze (`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`).

Freezes high-confidence clone definitions for primary analysis and
documents excluded/ambiguous calls, following the same mechanical,
SHA256-manifested freeze pattern already established in `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`
(`preprocess/release_freeze.py`) -- this step makes no new scientific
decisions of its own; Phases 8.00-8.07 already made every call (patient
mapping, false-positive estimation, resolution classification).

**High-confidence criterion, a judgment call made and documented
directly:** a clone (`08_tcr_clonal_analysis/07_build_clone_metadata_table.py`) is `high_confidence` if it
has at least one `singlet`-resolution cell (`08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`) -- i.e. is not
built entirely from `low_confidence` detections. Checked
against the data before finalising: this is equivalent here to "majority
of the clone's cells are singlet" (both criteria give the identical
166/213 split), so the simpler, more conservative "any singlet cell"
version is used without loss of strictness. The excluded 47 clones (built
entirely from `low_confidence` cells, i.e. every detection involved a
probe with `empirical_fpr > 0.5`, `08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`, `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`) are retained in the frozen
release's audit trail, not silently dropped from the record.

**Release status:** `CONDITIONAL GO`, matching the same framing already
used for the `04_quality_control/09_generate_qc_release_report.py` QC release -- multiple caveats accumulate
across TCR Clonal Analysis (only 105/216 probes have a statistically
identified intended patient, `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`; a tail of high-false-positive-rate
probes, `08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`; unresolved TRA/TRB pairing,
`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`, `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`, `08_tcr_clonal_analysis/07_build_clone_metadata_table.py`) that any downstream consumer of this release must
be aware of, not a clean unconditional release.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

RELEASE_NAME = "v1_tcr_calls"

SUPPORTING_FILES = [
    "metadata/tcr_probe_registry.tsv",
    "metadata/clone_ascertainment.tsv",
    "reports/tcr/patient_probe_audit.tsv",
    "data/derived/tcr_false_positive_estimates.parquet",
    "reports/tcr/cdr3_similarity_screen.tsv",
    "data/derived/tcr_resolved_calls.parquet",
]


def compute_file_hash(path: Path) -> str:
    sha256 = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def build_go_no_go_decision(clone_metadata: pd.DataFrame, ascertainment: pd.DataFrame) -> dict:
    fraction_identified = ascertainment["intended_patient_identified"].mean()
    fraction_high_confidence_clones = (clone_metadata["n_singlet_cells"] > 0).mean()
    caveats = [
        f"{(1 - fraction_identified) * 100:.1f}% of probes have no statistically identified intended patient (`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`).",
        "TRA/TRB pairing within a cell is inferred (`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s likely_single_clone_tra_trb_pair), not confirmed by any ground-truth linkage data (`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`).",
        f"{(1 - fraction_high_confidence_clones) * 100:.1f}% of identified clones are built entirely from low_confidence detections and excluded from the high-confidence set.",
        "Clone ascertainment (which clonotypes were probed at all) is abundance-biased with deliberate bystander/viral-reactive inclusions, not a random repertoire sample (`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`).",
    ]
    return {"status": "CONDITIONAL GO", "caveats": caveats}


def freeze_tcr_calls(project_root: Path, release_dir: Path) -> dict:
    clone_metadata_path = project_root / "data" / "derived" / "clone_metadata.parquet"
    ascertainment_path = project_root / "metadata" / "clone_ascertainment.tsv"

    if not clone_metadata_path.is_file():
        raise PipelineError(
            f"'{clone_metadata_path}' not found. Run `08_tcr_clonal_analysis/07_build_clone_metadata_table.py` first."
        )
    if not ascertainment_path.is_file():
        raise PipelineError(
            f"'{ascertainment_path}' not found. Run `08_tcr_clonal_analysis/02_document_clone_ascertainment.py` first."
        )

    clone_metadata = pd.read_parquet(clone_metadata_path)
    ascertainment = pd.read_csv(ascertainment_path, sep="\t")

    high_confidence = clone_metadata[clone_metadata["n_singlet_cells"] > 0].copy()
    excluded = clone_metadata[clone_metadata["n_singlet_cells"] == 0].copy()

    decision = build_go_no_go_decision(clone_metadata, ascertainment)

    release_dir.mkdir(parents=True, exist_ok=True)
    high_confidence_path = release_dir / "high_confidence_clones.parquet"
    excluded_path = release_dir / "excluded_low_confidence_clones.parquet"
    high_confidence.to_parquet(high_confidence_path)
    excluded.to_parquet(excluded_path)

    copied_files = [high_confidence_path, excluded_path]
    for rel_path in SUPPORTING_FILES:
        source = project_root / rel_path
        if source.is_file():
            dest = release_dir / Path(rel_path).name
            shutil.copy2(source, dest)
            copied_files.append(dest)

    checksums = {f.name: compute_file_hash(f) for f in copied_files}
    manifest = {
        "release_name": RELEASE_NAME,
        "frozen_at": datetime.now(timezone.utc).isoformat(),
        "status": decision["status"],
        "caveats": decision["caveats"],
        "n_high_confidence_clones": len(high_confidence),
        "n_excluded_clones": len(excluded),
        "n_high_confidence_cells": int(high_confidence["n_cells"].sum()),
        "files": checksums,
    }
    manifest_path = release_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False))

    checksums_path = release_dir / "checksums.sha256"
    checksums_path.write_text("".join(f"{h}  {name}\n" for name, h in checksums.items()))

    return {
        "release_dir": str(release_dir),
        "status": decision["status"],
        "n_high_confidence_clones": len(high_confidence),
        "n_excluded_clones": len(excluded),
        "n_high_confidence_cells": int(high_confidence["n_cells"].sum()),
        "n_files": len(copied_files),
    }
