#!/usr/bin/env python3
"""
`10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py`

Quantifies six per-section ecosystem properties (abundance,
fragmentation, interface length, mixing, compactness, relation to
tumour borders) from `10_niche_and_ecosystem_discovery/03_segment_tissue_domains.py`, `10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`'s domains and labels, Phase
9.03's calibrated graph, and Tumour Epithelium Characterisation's tumour boundary field -- see
src/xenium_tcr_ecology/niches/ecosystem_metrics.py's module docstring
.

Primary output: data/derived/ecosystem_metrics.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.niches.ecosystem_metrics import build_ecosystem_metrics


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
        script_name="05_quantify_ecosystem_abundance_and_topology",
        project_root=project_root,
        phase="10_niche_and_ecosystem_discovery",
    )

    try:
        summary = build_ecosystem_metrics(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_rows']} (section, ecosystem) row(s), {summary['n_sections']} section(s), "
        f"{summary['n_ecosystems']} ecosystem(s). Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
