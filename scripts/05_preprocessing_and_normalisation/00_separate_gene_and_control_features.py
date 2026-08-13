#!/usr/bin/env python3
"""
`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`

Builds a complete feature dictionary across all 18 sections' raw
transcripts tables (not just the 623 genes in the analysis matrix),
classifying every feature into an explicit class: biological_gene,
hpv_probe, tcr_cdr3_probe, negative_control_probe,
negative_control_codeword, or unassigned_codeword.

Primary output: metadata/feature_annotation.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.preprocess.feature_classification import build_feature_annotation_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    output_path = project_root / "metadata" / "feature_annotation.tsv"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "05_preprocessing_and_normalisation",
        script_name="00_separate_gene_and_control_features",
        project_root=project_root,
        phase="05_preprocessing_and_normalisation",
    )

    try:
        summary = build_feature_annotation_report(spatialdata_root, output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_features_total']} feature(s) across {summary['n_sections']} section(s). "
        f"Class counts: {summary['class_counts']}. "
        f"{summary['n_patient_specific_features']} patient-specific feature(s). "
        f"Wrote {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
