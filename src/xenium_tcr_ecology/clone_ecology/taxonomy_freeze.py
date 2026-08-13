"""Provisional taxonomy version freeze (`11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`).

Records `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`, `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s structure-test outcome as an
explicit, idempotent entry in `governance/taxonomy_version_log.tsv` --
`taxonomy_version = "v1_provisional"`. Provisional, not a final
commitment: External Checkpoint Validation ("Early external sanity-check and taxonomy-freeze
gate", `scripts/12_external_checkpoint_validation/`) projects this structure onto an
independent public reference and only then decides freeze-or-revise
(12.03/12.04); this milestone's `status` field is recorded as
`pending_phase12_external_validation`, not `frozen`.

Reads `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s decision (`clone_structure_test_results.
parquet`'s `final_decision`) rather than assuming a structure type --
this is currently `"continuous"`, so this module implements the continuous
branch, reading `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s factor-score output and loadings. An
explicit `PipelineError` guards against silently mishandling a future
re-run where `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` concludes `"discrete"` instead -- that
would need this module extended with a discrete-branch summary, not
silently forced through the continuous path.

Idempotent by `taxonomy_version`: re-running this script (e.g. after
upstream data changes) updates the existing `v1_provisional` row rather
than appending a duplicate -- "freeze" names a specific, current version
identity, not an ever-growing history log (that is what git history and
`results/logs/` are for).
"""

from __future__ import annotations

import datetime
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

TAXONOMY_VERSION = "v1_provisional"
FREEZE_STATUS = "pending_phase12_external_validation"


def build_taxonomy_version_entry(
    structure_type: str,
    n_units: int,
    dominant_feature: str,
    dominant_loading: float,
    frozen_date: str,
) -> dict:
    """Pure, testable: assembles one taxonomy-version-log row."""
    return {
        "taxonomy_version": TAXONOMY_VERSION,
        "structure_type": structure_type,
        "n_clone_sections": n_units,
        "dominant_feature": dominant_feature,
        "dominant_loading": dominant_loading,
        "source_script": "11_clone_spatial_descriptors/06_discover_provisional_structure.R",
        "status": FREEZE_STATUS,
        "frozen_date": frozen_date,
    }


def merge_taxonomy_log(existing: pd.DataFrame | None, new_entry: dict) -> pd.DataFrame:
    """Pure, testable: replaces any existing row with the same
    `taxonomy_version` (idempotent freeze) rather than duplicating it."""
    new_row = pd.DataFrame([new_entry])
    if existing is None or len(existing) == 0:
        return new_row
    kept = existing[existing["taxonomy_version"] != new_entry["taxonomy_version"]]
    return pd.concat([kept, new_row], ignore_index=True)


def build_taxonomy_freeze(project_root: Path) -> dict:
    structure_test_path = project_root / "data" / "derived" / "clone_structure_test_results.parquet"
    ecological_structure_path = (
        project_root / "data" / "derived" / "clone_ecological_structure.parquet"
    )
    loadings_path = (
        project_root / "data" / "derived" / "clone_ecological_structure_loadings.parquet"
    )
    output_path = project_root / "governance" / "taxonomy_version_log.tsv"

    for path, phase in [
        (
            structure_test_path,
            "`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`",
        ),
        (
            ecological_structure_path,
            "`11_clone_spatial_descriptors/06_discover_provisional_structure.R`",
        ),
        (loadings_path, "`11_clone_spatial_descriptors/06_discover_provisional_structure.R`"),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    structure_test = pd.read_parquet(structure_test_path)
    final_decisions = structure_test["final_decision"].unique()
    if len(final_decisions) != 1:
        raise PipelineError(
            f"`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s 'final_decision' is not uniform across analyses ({final_decisions.tolist()}) -- cannot freeze an ambiguous structure type."
        )
    structure_type = str(final_decisions[0])

    if structure_type != "continuous":
        raise PipelineError(
            f"`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` concluded structure_type='{structure_type}', but this module only implements the "
            "'continuous' branch (`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`'s current result). Extend build_taxonomy_freeze with "
            "a discrete-branch summary before running this on a discrete-structure result."
        )

    ecological_structure = pd.read_parquet(ecological_structure_path)
    loadings = pd.read_parquet(loadings_path)
    dominant = loadings.iloc[0]

    entry = build_taxonomy_version_entry(
        structure_type=structure_type,
        n_units=len(ecological_structure),
        dominant_feature=str(dominant["feature"]),
        dominant_loading=float(dominant["loading"]),
        frozen_date=datetime.date.today().isoformat(),
    )

    existing = pd.read_csv(output_path, sep="\t") if output_path.is_file() else None
    result = merge_taxonomy_log(existing, entry)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "taxonomy_version": TAXONOMY_VERSION,
        "structure_type": structure_type,
        "n_clone_sections": entry["n_clone_sections"],
        "dominant_feature": entry["dominant_feature"],
        "status": FREEZE_STATUS,
        "output_path": str(output_path),
    }
