"""Exports non-sensitive processed/derived data, with a licence and a
hash-manifested provenance record, for public archival (`17_statistical_
closure_and_release/06_build_public_data_release.py`).

Primary output: `release/data/`. Distinct from `release/software/`
(`07_build_documented_software_package.py`'s output) and from
`data/releases/final_primary/` (this script's own primary source,
produced by `00_freeze_primary_results.py`) -- this is the public-facing
bundle assembled *from* those internal releases, not a replacement for
either.

**Scope:** only aggregate, already-manuscript-facing outputs. This excludes anything
cell-level or per-patient beyond the de-identified patient/section codes
already public in the manuscript's own tables. `data/derived/
final_cell_annotations.parquet`, for example, is never included here,
regardless of its own privacy status, because it is not one of the three
source categories below:

1. `data/releases/final_primary/` -- the frozen, checksummed primary
   confirmatory results (Q1-Q3, HPV contrast), copied under
   `primary_results/` with its own nested MANIFEST.json/checksums.sha256
   preserved verbatim (that freeze's own self-contained provenance
   record, not superseded by this script's top-level one).
2. `tables/Table_*.tsv` -- the 10 manuscript result tables, copied under
   `tables/`. These are the same tables the manuscript and its
   Supplementary Information already cite; nothing here is disclosed
   for the first time by this export.
3. `metadata/data_dictionary.xlsx` -- the single schema reference for
   every derived field, copied under `metadata/`.

This script does not, by itself, satisfy the privacy/licensing review
the release procedure requires. It assembles exactly the
already-manuscript-facing scope above and nothing beyond it; that scope
is confirmed before this directory is uploaded anywhere.
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

PRIMARY_RESULTS_SOURCE_FILES = [
    "data/releases/final_primary/hpv_primary_contrasts.yaml",
    "data/releases/final_primary/hpv_robustness_summary.parquet",
    "data/releases/final_primary/clone_structure_test_results.parquet",
    "data/releases/final_primary/barrier_topology_model_results.parquet",
    "data/releases/final_primary/hpv_composition_comparison_results.parquet",
    "data/releases/final_primary/framework_generalisation_results.parquet",
    "data/releases/final_primary/variance_partition_results.parquet",
    "data/releases/final_primary/hpv_structure_comparison_results.parquet",
    "data/releases/final_primary/MANIFEST.json",
    "data/releases/final_primary/checksums.sha256",
]

TABLE_SOURCE_FILES = [
    "tables/Table_1_sample_manifest.tsv",
    "tables/Table_2_replicate_concordance.tsv",
    "tables/Table_3_statistical_summary.tsv",
    "tables/Table_4_variance_partition_results.tsv",
    "tables/Table_5_barrier_topology_model_results.tsv",
    "tables/Table_6_clone_structure_test_results.tsv",
    "tables/Table_7_claim_evidence_matrix.tsv",
    "tables/Table_8_source_paper_comparison.tsv",
    "tables/Table_9_validation_plan.tsv",
    "tables/Table_10_hpv_claim_strength.tsv",
]

METADATA_SOURCE_FILES = [
    "metadata/data_dictionary.xlsx",
]

# (source relative path, destination subdirectory under release/data/)
ALL_SOURCE_FILES = (
    [(f, "primary_results") for f in PRIMARY_RESULTS_SOURCE_FILES]
    + [(f, "tables") for f in TABLE_SOURCE_FILES]
    + [(f, "metadata") for f in METADATA_SOURCE_FILES]
)

LICENSE_TEXT = """This data release is licensed under the Creative Commons
Attribution 4.0 International License (CC BY 4.0).

You are free to share and adapt this material for any purpose, including
commercially, as long as you give appropriate credit, provide a link to
the license, and indicate if changes were made.

Full legal text: https://creativecommons.org/licenses/by/4.0/legalcode

Cite the software DOI (see this repository's CITATION.cff) and the
manuscript this data accompanies.
"""


def _build_readme(manifest_analysis_ids: list[str]) -> str:
    return f"""# xenium-tcr-ecology: public data release

Processed, non-sensitive derived data underlying the accompanying
manuscript's reported results. See `MANIFEST.json` for the exact file
list and SHA-256 checksums, and `LICENSE` for reuse terms.

## Contents

- `primary_results/` -- the frozen, prespecified primary confirmatory
  results ({", ".join(manifest_analysis_ids)}), with the original
  freeze's own nested `MANIFEST.json`/`checksums.sha256` preserved
  verbatim (see `xenium_tcr_ecology.release.freeze_primary_results` in
  the source code archive for how that freeze was produced).
- `tables/` -- the 10 manuscript result tables (`Table_1`..`Table_10`),
  the same tables cited throughout the manuscript and its Supplementary
  Information.
- `metadata/data_dictionary.xlsx` -- the schema reference for every
  derived field appearing in the files above.

## What this release excludes

No cell-level or per-cell spatial data, and no data beyond the
de-identified patient/section codes (e.g. `P01`, `P01_run1`) already
public in the manuscript's own tables, is included here. Raw and
per-cell derived data are not redistributed by this release; the raw
sequencing/imaging data are available from their original public
accessions (GEO, 10x Genomics, UCSC Xena -- see the manuscript's
"Availability of data and materials" and `manifests/dataset_registry.yaml`
in the source code archive), and the full computational workflow that
produces every intermediate and per-cell artefact from that raw data is
in the accompanying source code release (see the manuscript's
"Availability of data and materials" section, or this deposit's own
Zenodo metadata, for the exact software version/DOI this data was built
from).

## Licence

CC BY 4.0 -- see `LICENSE`.

## Provenance

Generated by `17_statistical_closure_and_release/
06_build_public_data_release.py` in the accompanying source code
repository. Re-running that script against the same frozen inputs
reproduces this directory's contents and checksums exactly.
"""


def build_public_data_release(project_root: Path, release_dir: Path | None = None) -> dict:
    release_dir = release_dir or (project_root / "release" / "data")

    missing = [f for f, _ in ALL_SOURCE_FILES if not (project_root / f).is_file()]
    if missing:
        raise PipelineError(
            f"Missing public-data-release source file(s), cannot build release/data/: {missing}. "
            "Run `00_freeze_primary_results.py`, `04_generate_results_tables.py`, and "
            "`02_generate_data_dictionary.py` first."
        )

    current_checksums = {
        f"{subdir}/{Path(f).name}": compute_file_hash(project_root / f)
        for f, subdir in ALL_SOURCE_FILES
    }

    existing_manifest_path = release_dir / "MANIFEST.json"
    if existing_manifest_path.is_file():
        previous_manifest = json.loads(existing_manifest_path.read_text())
        previous_checksums = previous_manifest.get("files", {})
        changed = check_hash_consistency(previous_checksums, current_checksums)
        if changed:
            raise PipelineError(
                f"Upstream input(s) changed since 'release/data/' was last built: {changed}. "
                "Refusing to silently rebuild the public data release over changed content -- "
                "confirm the change is intended (e.g. a corrected result) before rebuilding, "
                "and re-review privacy/licensing before re-publishing."
            )

    release_dir.mkdir(parents=True, exist_ok=True)
    copied_files: list[Path] = []
    for rel_path, subdir in ALL_SOURCE_FILES:
        source = project_root / rel_path
        dest_dir = release_dir / subdir
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / Path(rel_path).name
        dest.write_bytes(source.read_bytes())
        copied_files.append(dest)

    checksums = {str(f.relative_to(release_dir)): compute_file_hash(f) for f in copied_files}

    primary_manifest_path = project_root / "data" / "releases" / "final_primary" / "MANIFEST.json"
    primary_analysis_ids = json.loads(primary_manifest_path.read_text())["primary_analysis_ids"]

    manifest = {
        "release_name": "public_data_release",
        "built_at": datetime.now(timezone.utc).isoformat(),
        "license": "CC-BY-4.0",
        "primary_analysis_ids": primary_analysis_ids,
        "excludes": (
            "cell-level and per-cell spatial data; anything beyond de-identified "
            "patient/section codes already public in the manuscript's own tables"
        ),
        "files": checksums,
    }
    manifest_path = release_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False))

    checksums_path = release_dir / "checksums.sha256"
    checksums_path.write_text("".join(f"{h}  {name}\n" for name, h in checksums.items()))

    (release_dir / "LICENSE").write_text(LICENSE_TEXT)
    (release_dir / "README.md").write_text(_build_readme(primary_analysis_ids))

    return {
        "release_dir": str(release_dir),
        "n_files": len(copied_files),
        "license": "CC-BY-4.0",
    }
