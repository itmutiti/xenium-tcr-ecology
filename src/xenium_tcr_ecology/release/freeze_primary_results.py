"""Freezes the prespecified primary results (Q1-Q3, plus the single
HPV contrast set) into an immutable, hash-manifested release directory
before any exploratory extension work builds on them (`17_statistical_closure_and_release/00_freeze_primary_results.py`).

**Direct reuse of the established freeze pattern, not a fresh design:**
`compute_file_hash` and `check_hash_consistency` are imported directly
from `13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py`'s `clone_ecology.taxonomy_release` (itself following
the same mechanical pattern already established in `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`, `08_tcr_clonal_analysis/08_generate_tcr_release_report.py`) --
this milestone makes no new scientific decisions, it locks in decisions
Phases 9, 11-16 already made. On the first run there is no prior
release to compare against, so this establishes the reference hashes;
on every subsequent run, if any upstream input has changed, this raises
a `PipelineError` rather than silently re-freezing over a changed
scientific result.

**Exact scope: only the primary family, not every registered
analysis.** `governance/analysis_registry.tsv`'s `multiplicity_
family` column marks exactly 5 entries `primary`:
`q1_framework_generalisation`, `q2_variance_partition_confirmatory`,
`q2_discrete_vs_continuous_structure_test`, `q3_barrier_topology_
confirmatory`, `hpv_primary_contrast`. The `supporting`/`methods
validation`/robustness-check entries (`q3_literature_benchmark`,
`null_model_calibration_suite`, `tcr_false_positive_rate`,
`replicate_concordance`, `segmentation_robustness_check`, `external_
checkpoint_directional_consistency`) are important evidence (`16_external_validation_and_generalisation/07_generate_evidence_matrix.py`'s evidence matrix links every one of them) but are not
part of the prespecified primary family this milestone's scaffold
names, and are deliberately excluded from this specific freeze.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


from xenium_tcr_ecology.clone_ecology.taxonomy_release import (
    check_hash_consistency,
    compute_file_hash,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError

RELEASE_NAME = "final_primary"

PRIMARY_RESULT_FILES = [
    # q1_framework_generalisation
    "data/derived/framework_generalisation_results.parquet",
    # q2_variance_partition_confirmatory
    "data/derived/variance_partition_results.parquet",
    # q2_discrete_vs_continuous_structure_test
    "data/derived/clone_structure_test_results.parquet",
    # q3_barrier_topology_confirmatory
    "data/derived/barrier_topology_model_results.parquet",
    # hpv_primary_contrast
    "governance/hpv_primary_contrasts.yaml",
    "data/derived/hpv_composition_comparison_results.parquet",
    "data/derived/hpv_structure_comparison_results.parquet",
    "data/derived/hpv_robustness_summary.parquet",
]


def build_primary_results_freeze(project_root: Path, release_dir: Path | None = None) -> dict:
    release_dir = release_dir or (project_root / "data" / "releases" / RELEASE_NAME)

    missing = [f for f in PRIMARY_RESULT_FILES if not (project_root / f).is_file()]
    if missing:
        raise PipelineError(
            f"Missing prespecified primary-result input file(s), cannot freeze: {missing}"
        )

    current_checksums = {
        Path(f).name: compute_file_hash(project_root / f) for f in PRIMARY_RESULT_FILES
    }

    existing_manifest_path = release_dir / "MANIFEST.json"
    if existing_manifest_path.is_file():
        previous_manifest = json.loads(existing_manifest_path.read_text())
        previous_checksums = previous_manifest.get("files", {})
        changed = check_hash_consistency(previous_checksums, current_checksums)
        if changed:
            raise PipelineError(
                f"Upstream input(s) changed since '{RELEASE_NAME}' was last frozen: {changed}. "
                "Refusing to silently re-freeze over a changed scientific result."
            )

    release_dir.mkdir(parents=True, exist_ok=True)
    copied_files = []
    for rel_path in PRIMARY_RESULT_FILES:
        source = project_root / rel_path
        dest = release_dir / Path(rel_path).name
        dest.write_bytes(source.read_bytes())
        copied_files.append(dest)

    checksums = {f.name: compute_file_hash(f) for f in copied_files}

    manifest = {
        "release_name": RELEASE_NAME,
        "frozen_at": datetime.now(timezone.utc).isoformat(),
        "primary_analysis_ids": [
            "q1_framework_generalisation",
            "q2_variance_partition_confirmatory",
            "q2_discrete_vs_continuous_structure_test",
            "q3_barrier_topology_confirmatory",
            "hpv_primary_contrast",
        ],
        "files": checksums,
    }
    manifest_path = release_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False))

    checksums_path = release_dir / "checksums.sha256"
    checksums_path.write_text("".join(f"{h}  {name}\n" for name, h in checksums.items()))

    return {
        "release_dir": str(release_dir),
        "n_files": len(copied_files),
        "n_primary_analyses": len(manifest["primary_analysis_ids"]),
    }
