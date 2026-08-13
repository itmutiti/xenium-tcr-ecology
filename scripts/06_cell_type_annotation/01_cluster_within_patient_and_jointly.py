#!/usr/bin/env python3
"""
`06_cell_type_annotation/01_cluster_within_patient_and_jointly.py`

Generates Leiden clustering at multiple resolutions on the pooled dataset
("jointly") plus a single-resolution clustering computed independently
within each patient, restricted to biological_gene features. Cluster
labels are exploratory structure, not cell-type calls -- see
src/xenium_tcr_ecology/annotation/clustering.py's module docstring.

Primary output: data/derived/clustering_assignments.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.clustering import build_clustering_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "06_cell_type_annotation",
        script_name="01_cluster_within_patient_and_jointly",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_clustering_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells, {summary['n_clustering_genes']} clustering genes. "
        f"Clusters by resolution: {summary['n_clusters_by_resolution']}. "
        f"{summary['n_within_patient_clusters']} within-patient cluster(s) across "
        f"{summary['n_cells_with_within_patient_label']:,} cells. "
        f"Wrote data/derived/clustering_assignments.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
