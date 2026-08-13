#!/usr/bin/env python3
"""
`10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`

Assigns a descriptive ecosystem name to each of `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s
unsupervised archetypes via a documented, mechanical enrichment-ratio
rubric applied AFTER discovery -- see
src/xenium_tcr_ecology/niches/ecosystem_annotation.py's module docstring

this differs from `06_cell_type_annotation/07_blinded_annotation_review.py`'s blinded pathologist review.

Primary output: metadata/ecosystem_annotation.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.niches.ecosystem_annotation import build_ecosystem_annotation


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "10_niche_and_ecosystem_discovery",
        script_name="04_annotate_ecosystems_with_blinded_rules",
        project_root=project_root,
        phase="10_niche_and_ecosystem_discovery",
    )

    try:
        summary = build_ecosystem_annotation(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_archetypes']} archetype(s), {summary['n_mixed_non_specific']} mixed/non-specific. "
        f"Labels: {summary['ecosystem_labels']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
