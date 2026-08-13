#!/usr/bin/env python3
"""
`06_cell_type_annotation/03_map_external_scrna_reference.py`

Maps GSE287301 (McCord et al. 2026's own companion scRNA-seq dataset --
T cells from the same 28-patient HNSCC cohort) onto the Xenium data via a
nearest-centroid Pearson-correlation classifier restricted to
Xenium-panel-overlapping genes. Reference T-cell-state labels are derived
independently via marker-based scoring (the paper's own cluster labels are
locked in a proprietary Loupe Browser file, not separately available as a
scriptable table). Transfer-confidence degradation from the restricted gene
overlap is explicitly benchmarked, not asserted -- see
src/xenium_tcr_ecology/annotation/reference_mapping.py's module docstring.

Primary output: data/derived/reference_labels.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.reference_mapping import build_reference_mapping_report


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
        script_name="03_map_external_scrna_reference",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_reference_mapping_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    deg = summary["transfer_confidence_degradation"]
    print(
        f"[OK]   Reference: {summary['reference_n_cells_after_qc']:,} cells, "
        f"{summary['n_t_cell_states']} T-cell states. "
        f"Degradation benchmark: broad-gene-set accuracy={deg['broad_gene_set_accuracy']:.3f}, "
        f"panel-restricted accuracy={deg['panel_restricted_accuracy']:.3f} "
        f"({deg['n_panel_overlap_genes']} overlap genes). "
        f"Xenium: {summary['n_xenium_cells_labeled']:,} cells labeled, "
        f"mean confidence={summary['xenium_mean_confidence']:.3f}. "
        f"Wrote data/derived/reference_labels.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
