#!/usr/bin/env python3
"""
`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`

Predeclares 4 biologically motivated sender-receiver comparisons
(tumour/fibroblast/myeloid/APC -> T cell), each grounded directly in
this project's already-established findings -- see
src/xenium_tcr_ecology/interactions/sender_receiver_pairs.py's module
docstring (including the zero-coverage TGF-beta panel gap).

Primary output: config/sender_receiver_pairs.yaml
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.interactions.sender_receiver_pairs import build_sender_receiver_config


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
        script_name="00_define_sender_receiver_pairs",
        project_root=project_root,
        phase="14_spatial_interactions_and_barriers",
    )

    try:
        summary = build_sender_receiver_config(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_pairs']} sender-receiver pair(s), {summary['n_programs']} programme gene set(s). "
        f"Zero panel coverage: {summary['programs_with_zero_coverage']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
