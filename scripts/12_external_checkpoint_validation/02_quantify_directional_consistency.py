#!/usr/bin/env python3
"""
`12_external_checkpoint_validation/02_quantify_directional_consistency.py`

Quantifies pairwise sign agreement and rank correlation of T-cell-
state abundances between this project and the independent GSE103322
reference, and links each state's rank shift to `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s
continuous-structure loadings -- see
src/xenium_tcr_ecology/external_checkpoint/directional_consistency.py's module
docstring.

Primary output: results/external_checkpoint/directional_consistency.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.external_checkpoint.directional_consistency import (
    build_directional_consistency,
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
        logs_dir=project_root / "results" / "logs" / "12_external_checkpoint_validation",
        script_name="02_quantify_directional_consistency",
        project_root=project_root,
        phase="12_external_checkpoint_validation",
    )

    try:
        summary = build_directional_consistency(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_pairs_agreeing']}/{summary['n_state_pairs']} state pairs directionally agree "
        f"({summary['fraction_pairs_agreeing']:.1%}). Spearman rank correlation rho={summary['spearman_rho']:.3f} "
        f"(p={summary['spearman_pvalue']:.3f}). Wrote {summary['output_path']}, {summary['sign_agreement_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
