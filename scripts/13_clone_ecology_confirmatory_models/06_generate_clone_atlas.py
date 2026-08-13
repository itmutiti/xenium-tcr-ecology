#!/usr/bin/env python3
"""
`13_clone_ecology_confirmatory_models/06_generate_clone_atlas.py`

Produces one spatial thumbnail plus a standardised descriptor summary
for every high-confidence clone-section, in a single,
systematically-ordered HTML page, preventing example-selection bias by
construction -- see
src/xenium_tcr_ecology/clone_ecology/clone_atlas.py's module docstring
.

Primary output: reports/clone_ecology/clone_atlas.html
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.clone_atlas import build_clone_atlas


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "13_clone_ecology_confirmatory_models",
        script_name="06_generate_clone_atlas",
        project_root=project_root,
        phase="13_clone_ecology_confirmatory_models",
    )

    try:
        summary = build_clone_atlas(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_clone_section_rows']} clone-section row(s), {summary['n_distinct_clones']} distinct clone(s), "
        f"{summary['n_sections']} section(s). Wrote {summary['output_path']}, {summary['merged_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
