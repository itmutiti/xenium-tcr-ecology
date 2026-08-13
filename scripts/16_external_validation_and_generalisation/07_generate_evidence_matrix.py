#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/07_generate_evidence_matrix.py`

Links every one of this project's registered claims (`governance/
analysis_registry.tsv`) to its discovery, sensitivity, replicate
and external-validation evidence, with an overall grade.
See src/xenium_tcr_ecology/validation/evidence_matrix.py's module
docstring.

Primary output: results/claim_evidence_matrix.tsv
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.evidence_matrix import build_claim_evidence_matrix


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "16_external_validation_and_generalisation",
        script_name="07_generate_evidence_matrix",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_claim_evidence_matrix(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_claims']} registered claim(s) linked to evidence: {summary['grade_counts']}. Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
