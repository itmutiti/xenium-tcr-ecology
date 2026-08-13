"""Frozen taxonomy version release + hash-change gate (`13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py`).

Freezes `taxonomy_version = v1_provisional` (`11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`, `12_external_checkpoint_validation/03_decide_freeze_or_revise.py`-12.04)
into a SHA256-manifested release directory, following the same
mechanical freeze pattern already established in `08_tcr_clonal_analysis/08_generate_tcr_release_report.py`, `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`
(`tcr/release_freeze.py`, `preprocess/release_freeze.py`) -- this step
makes no new scientific decisions; Clone Spatial Descriptors/12 already made every call.

**"Refuses to run if any upstream input hash has changed" -- idempotent
design:** on the first run, there is no prior release to compare
against, so this establishes the reference hashes (recorded in
`MANIFEST.json`, matching every other release in this project). On
every subsequent run, this recomputes the same files' hashes and
compares them against the existing release's recorded hashes -- if
`data/derived/clone_ecological_structure.parquet` (or any other
upstream input) has changed since the release was frozen, this raises a
`PipelineError` rather than silently re-freezing over a changed
scientific result. `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`-13.06's confirmatory models are meant to
run against one frozen, stable clone-structure definition -- allowing it
to silently drift between runs would undermine the entire point of a
freeze gate.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

RELEASE_NAME = "v1_clone_structure"

UPSTREAM_FILES = [
    "data/derived/clone_ecological_structure.parquet",
    "data/derived/clone_ecological_structure_loadings.parquet",
    "data/derived/clone_structure_test_results.parquet",
    "data/derived/clone_spatial_descriptors.parquet",
    "data/derived/clone_state_composition.parquet",
    "data/derived/clone_tumour_engagement.parquet",
    "data/derived/clone_apc_support.parquet",
    "data/derived/clone_barrier_metrics.parquet",
    "governance/taxonomy_version_log.tsv",
    "governance/freeze_decision.tsv",
]


def compute_file_hash(path: Path) -> str:
    sha256 = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def check_hash_consistency(
    previous_checksums: dict[str, str], current_checksums: dict[str, str]
) -> list[str]:
    """Pure, testable: file names whose hash differs between two
    checksum dicts (only compares files present in both -- a file added
    or removed between releases is not itself a "changed" file)."""
    changed = []
    for name, prev_hash in previous_checksums.items():
        current_hash = current_checksums.get(name)
        if current_hash is not None and current_hash != prev_hash:
            changed.append(name)
    return sorted(changed)


def build_taxonomy_release(project_root: Path, release_dir: Path | None = None) -> dict:
    release_dir = release_dir or (project_root / "data" / "releases" / RELEASE_NAME)

    taxonomy_log_path = project_root / "governance" / "taxonomy_version_log.tsv"
    freeze_decision_path = project_root / "governance" / "freeze_decision.tsv"
    if not taxonomy_log_path.is_file():
        raise PipelineError(
            f"'{taxonomy_log_path}' not found. Run `11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py` first."
        )
    if not freeze_decision_path.is_file():
        raise PipelineError(
            f"'{freeze_decision_path}' not found. Run `12_external_checkpoint_validation/03_decide_freeze_or_revise.py` first."
        )

    freeze_decision = pd.read_csv(freeze_decision_path, sep="\t").iloc[-1]
    if freeze_decision["decision"] != "freeze":
        raise PipelineError(
            f"`12_external_checkpoint_validation/03_decide_freeze_or_revise.py`'s decision was '{freeze_decision['decision']}', not 'freeze' -- "
            "this script only implements the freeze branch. A 'revise' decision requires a "
            "v2 taxonomy to be produced upstream before this release can be built."
        )

    missing = [f for f in UPSTREAM_FILES if not (project_root / f).is_file()]
    if missing:
        raise PipelineError(f"Missing upstream input file(s), cannot freeze: {missing}")

    current_checksums = {Path(f).name: compute_file_hash(project_root / f) for f in UPSTREAM_FILES}

    existing_manifest_path = release_dir / "MANIFEST.json"
    if existing_manifest_path.is_file():
        previous_manifest = json.loads(existing_manifest_path.read_text())
        previous_checksums = previous_manifest.get("files", {})
        changed = check_hash_consistency(previous_checksums, current_checksums)
        if changed:
            raise PipelineError(
                f"Upstream input(s) changed since '{RELEASE_NAME}' was last frozen: {changed}. "
                "Refusing to silently re-freeze over a changed scientific result -- investigate why "
                "these files changed before re-running this release."
            )

    release_dir.mkdir(parents=True, exist_ok=True)
    copied_files = []
    for rel_path in UPSTREAM_FILES:
        source = project_root / rel_path
        dest = release_dir / Path(rel_path).name
        shutil.copy2(source, dest)
        copied_files.append(dest)

    checksums = {f.name: compute_file_hash(f) for f in copied_files}
    structure = pd.read_parquet(release_dir / "clone_ecological_structure.parquet")

    manifest = {
        "release_name": RELEASE_NAME,
        "frozen_at": datetime.now(timezone.utc).isoformat(),
        "taxonomy_version": str(freeze_decision["taxonomy_version_evaluated"]),
        "decision": str(freeze_decision["decision"]),
        "flagged_states": str(freeze_decision["flagged_states"]),
        "n_clone_section_rows": len(structure),
        "n_distinct_clones": int(structure["clone_id"].nunique()),
        "files": checksums,
    }
    manifest_path = release_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False))

    checksums_path = release_dir / "checksums.sha256"
    checksums_path.write_text("".join(f"{h}  {name}\n" for name, h in checksums.items()))

    return {
        "release_dir": str(release_dir),
        "taxonomy_version": manifest["taxonomy_version"],
        "n_clone_section_rows": manifest["n_clone_section_rows"],
        "n_distinct_clones": manifest["n_distinct_clones"],
        "n_files": len(copied_files),
    }
