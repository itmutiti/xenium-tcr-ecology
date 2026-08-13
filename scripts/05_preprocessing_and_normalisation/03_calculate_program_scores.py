#!/usr/bin/env python3
"""
`05_preprocessing_and_normalisation/03_calculate_program_scores.py`

Computes curated cytotoxicity, exhaustion, activation, interferon,
proliferation, stress, EMT and antigen-presentation scores using scanpy's
score_genes on `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`'s primary normalisation layer. Several programs
have documented panel-coverage gaps relative to canonical literature
signatures -- see src/xenium_tcr_ecology/preprocess/program_scores.py's
module docstring.

Primary output: data/derived/program_scores.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.preprocess.program_scores import build_program_scores_report


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "05_preprocessing_and_normalisation",
        script_name="03_calculate_program_scores",
        project_root=project_root,
        phase="05_preprocessing_and_normalisation",
    )

    try:
        summary = build_program_scores_report(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_cells']:,} cells, {summary['n_programs']} program(s), "
        f"layer='{summary['layer_used']}'. Genes per program: {summary['genes_per_program']}. "
        f"Thin-coverage programs (documented panel gaps): {summary['thin_coverage_programs']}. "
        f"Wrote data/derived/program_scores.parquet"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
