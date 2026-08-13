#!/usr/bin/env python3
"""
`11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`

Records `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R`, `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s structure-test outcome as a
provisional `taxonomy_version = "v1_provisional"` entry, routing to
External Checkpoint Validation for external sanity-checking before any final freeze -- see
src/xenium_tcr_ecology/clone_ecology/taxonomy_freeze.py's module
docstring.

Primary output: governance/taxonomy_version_log.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.clone_ecology.taxonomy_freeze import build_taxonomy_freeze


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
        script_name="07_freeze_provisional_taxonomy_version",
        project_root=project_root,
        phase="11_clone_spatial_descriptors",
    )

    try:
        summary = build_taxonomy_freeze(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   Froze {summary['taxonomy_version']} ({summary['structure_type']}, "
        f"{summary['n_clone_sections']} clone-sections, dominant feature '{summary['dominant_feature']}'). "
        f"Status: {summary['status']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
