#!/usr/bin/env python3
"""One-time data migration: corrects the dtype of three obs columns
(`included_in_primary_hnscc_cohort`, `is_technical_replicate`,
`hpv_p16_positive`) in already-materialised h5ad files that were written
before src/xenium_tcr_ecology/io/metadata_join.py's boolean-casting fix.

The underlying values have always been correct -- this only corrects their
storage dtype (category of the literal strings "True"/"False" -> native
bool), so a full re-run of `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` onward is not needed; this migration
is a cheap, low-risk, one-time correction of already-materialised files,
not a scientific recomputation. Run once, not part of the regular pipeline
DAG (no Snakemake rule) -- re-running scripts/03_spatialdata_import/04_attach_metadata.py
from scratch would already produce correctly-typed output directly.

Usage: python tools/migrate_boolean_metadata_dtype.py --project-root <path>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.metadata_join import BOOLEAN_METADATA_FIELDS

TARGET_FILES = [
    "data/objects/hnscc_xenium_combined.h5ad",
    "data/objects/qc_filtered.h5ad",
    "data/objects/analysis_ready.h5ad",
]


def migrate_file(path: Path) -> dict:
    adata = ad.read_h5ad(path)
    fixed_fields = []
    for field in BOOLEAN_METADATA_FIELDS:
        if field not in adata.obs.columns:
            continue
        if adata.obs[field].dtype == bool:
            continue
        adata.obs[field] = adata.obs[field].astype(str) == "True"
        fixed_fields.append(field)
    if fixed_fields:
        adata.write_h5ad(path)
    return {"path": str(path), "fixed_fields": fixed_fields}


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--project-root", default=None)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    for rel_path in TARGET_FILES:
        path = project_root / rel_path
        if not path.is_file():
            print(f"[SKIP]  '{path}' does not exist.")
            continue
        result = migrate_file(path)
        if result["fixed_fields"]:
            print(f"[OK]    {path}: fixed {result['fixed_fields']}")
        else:
            print(f"[OK]    {path}: already correctly typed, no change.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
