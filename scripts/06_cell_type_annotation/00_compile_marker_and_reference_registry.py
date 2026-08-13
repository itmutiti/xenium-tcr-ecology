#!/usr/bin/env python3
"""
`06_cell_type_annotation/00_compile_marker_and_reference_registry.py`

Builds a versioned registry mapping the 399-gene biological_gene panel
to feasible cell identities, at the resolution the panel actually supports.
Confidence tiers reflect genuine lineage specificity in HNSCC tissue, not
merely gene presence -- see
src/xenium_tcr_ecology/annotation/marker_registry.py's module docstring for
the evidence (a complete pancreatic-islet-hormone panel, renal/hepatic/
melanocyte markers) that this is a generic commercial multi-tissue panel,
not an HNSCC-bespoke design.

Primary output: references/cell_type_marker_registry.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.annotation.marker_registry import build_marker_registry_report


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
        script_name="00_compile_marker_and_reference_registry",
        project_root=project_root,
        phase="06_cell_type_annotation",
    )

    try:
        summary = build_marker_registry_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Registry {summary['registry_version']}: {summary['n_identities']} identit(y/ies) "
        f"({summary['n_major_lineages']} major lineage(s), {summary['n_substates']} substate(s)). "
        f"Confidence tiers: {summary['confidence_tier_counts']}. "
        f"{summary['n_panel_genes_mapped']}/{summary['n_panel_genes_total']} panel genes mapped "
        f"({summary['n_panel_genes_unmapped']} unmapped). "
        f"Wrote references/cell_type_marker_registry.tsv"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
