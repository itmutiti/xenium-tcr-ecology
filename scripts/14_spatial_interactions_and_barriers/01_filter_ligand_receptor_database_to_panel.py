#!/usr/bin/env python3
"""
`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`

Checks a hand-curated set of canonical immune checkpoint/
costimulatory/homing ligand-receptor pairs against the panel,
labelling incomplete pairs (only one side present) to prevent
overinterpretation -- see
src/xenium_tcr_ecology/interactions/ligand_receptor_database.py's module
docstring (including why a full generic LR database was not used for
this targeted panel).

Primary output: references/panel_supported_interactions.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.interactions.ligand_receptor_database import (
    build_panel_supported_interactions,
)


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "14_spatial_interactions_and_barriers",
        script_name="01_filter_ligand_receptor_database_to_panel",
        project_root=project_root,
        phase="14_spatial_interactions_and_barriers",
    )

    try:
        summary = build_panel_supported_interactions(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_pairs_total']} candidate pair(s): {summary['n_pairs_complete']} panel-complete, "
        f"{summary['n_pairs_incomplete']} incomplete. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
