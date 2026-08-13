#!/usr/bin/env python3
"""
`13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py`

Freezes taxonomy_version=v1_provisional (`11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`, `12_external_checkpoint_validation/03_decide_freeze_or_revise.py`) into a
SHA256-manifested release directory; refuses to re-freeze if any
upstream input has changed since the last freeze -- see
src/xenium_tcr_ecology/clone_ecology/taxonomy_release.py's module
docstring.

Primary output: data/releases/v1_clone_structure/
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.taxonomy_release import build_taxonomy_release


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
        script_name="00_load_frozen_taxonomy_version",
        project_root=project_root,
        phase="13_clone_ecology_confirmatory_models",
    )

    try:
        summary = build_taxonomy_release(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Froze {summary['taxonomy_version']}: {summary['n_clone_section_rows']} clone-section rows, "
        f"{summary['n_distinct_clones']} distinct clones, {summary['n_files']} files. Wrote {summary['release_dir']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
