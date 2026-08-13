#!/usr/bin/env python3
"""
`11_clone_spatial_descriptors/01_compute_clone_cell_state_composition.py`

Summarises each clone-section's T-cell state composition (Phase
6.04's taxonomy: Cytotoxic, Exhausted, Cycling, Treg, CD4, CD8,
Ambiguous), each with a Clopper-Pearson exact binomial confidence
interval, plus overall Shannon-entropy heterogeneity -- see
src/xenium_tcr_ecology/clone_ecology/state_composition.py's module
docstring.

Primary output: data/derived/clone_state_composition.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.state_composition import build_clone_state_composition


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "11_clone_spatial_descriptors",
        script_name="01_compute_clone_cell_state_composition",
        project_root=project_root,
        phase="11_clone_spatial_descriptors",
    )

    try:
        summary = build_clone_state_composition(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clone_section_rows']} (clone, section) row(s), "
        f"{summary['n_distinct_clones']} distinct clone(s). Mean Shannon entropy "
        f"{summary['mean_shannon_entropy_bits']:.3f} bits. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
