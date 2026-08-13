#!/usr/bin/env python3
"""
`03_spatialdata_import/04_attach_clinical_and_technical_metadata.py`

Joins patient/run metadata from metadata/sample_manifest.tsv onto each
section's AnnData (in place, overwriting the `03_spatialdata_import/03_create_anndata_expression_objects.py` output), validating
the join key (section_id) is unique on both sides before merging.

Primary output: updated h5ad objects; join audit (this script's log)
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.metadata_join import attach_metadata


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    anndata_root = project_root / "data" / "objects" / "anndata"
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    h5ad_paths = sorted(anndata_root.glob("*.h5ad"))
    if not h5ad_paths:
        print(
            f"[ERROR] No .h5ad files found under '{anndata_root}'. Run `03_spatialdata_import/03_create_anndata_expression_objects.py` first.",
            file=sys.stderr,
        )
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="04_attach_clinical_and_technical_metadata",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    joined, failed = 0, []
    for h5ad_path in h5ad_paths:
        try:
            summary = attach_metadata(h5ad_path, sample_manifest_path, h5ad_path)
        except PipelineError as exc:
            print(f"[ERROR] {h5ad_path.stem}: {exc}", file=sys.stderr)
            logger.log_error(f"{h5ad_path.stem}: {exc}")
            failed.append(h5ad_path.stem)
            continue
        print(
            f"[OK]   {summary['section_id']} (patient {summary['patient_id']}): {summary['n_cells']} cells"
        )
        logger.log_event(**summary)
        joined += 1

    status = "ok" if not failed else "failed"
    logger.write(status=status)
    print(f"\n[SUMMARY] {joined} joined, {len(failed)} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
