#!/usr/bin/env python3
"""
`01_project_setup_and_governance/01_build_sample_manifest.py`

Compiles config/geo/sample_manifest_input.yaml (per-sample metadata fetched
directly from GSE300147's 18 individual GSM SOFT-format records on GEO) into
the canonical metadata/sample_manifest.tsv: patient/section identifiers,
technical-replicate and primary-cohort flags, p16/HPV status, and clinical
covariates.

Scope note: file paths and checksums are NOT included here, despite being
named in this script's original description -- they do not exist until
Raw Data Ingestion downloads and verifies the raw archive. `02_raw_data_ingestion/00_query_geo_accession.py`-2.02's outputs
are joined onto this manifest (by gsm_accession) once files are
on disk, rather than this script fabricating placeholder paths now.

Primary output: metadata/sample_manifest.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.ingest.manifest import compile_sample_manifest


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    input_path = project_root / "config" / "geo" / "sample_manifest_input.yaml"
    output_path = project_root / "metadata" / "sample_manifest.tsv"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "01_project_setup_and_governance",
        script_name="01_build_sample_manifest",
        project_root=project_root,
        phase="01_project_setup_and_governance",
    )

    try:
        summary = compile_sample_manifest(input_path, output_path, project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(f"[OK]   Wrote {summary['total_samples']} sample(s) to {output_path}")
    print(
        f"[INFO] {summary['hnscc_patients']} HNSCC patients "
        f"({summary['replicated_patients']} with technical replicates, "
        f"{summary['hnscc_sections']} sections total), "
        f"{summary['hpv_positive_patients']} p16/HPV-positive, "
        f"{summary['ameloblastoma_specimens']} non-HNSCC specimen(s) excluded from primary cohort."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
