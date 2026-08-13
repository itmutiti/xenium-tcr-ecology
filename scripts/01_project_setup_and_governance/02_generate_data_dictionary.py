#!/usr/bin/env python3
"""
`01_project_setup_and_governance/02_generate_data_dictionary.py`

Compiles config/metadata/data_dictionary_input.yaml into
metadata/data_dictionary.xlsx (one worksheet per documented table), and
cross-checks each documented table against its on-disk columns where
that table already exists, so the dictionary cannot silently drift out of
sync with what earlier scripts actually produced.

Primary output: metadata/data_dictionary.xlsx
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.metadata.data_dictionary import compile_data_dictionary


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    input_path = project_root / "config" / "metadata" / "data_dictionary_input.yaml"
    output_path = project_root / "metadata" / "data_dictionary.xlsx"
    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "01_project_setup_and_governance",
        script_name="02_generate_data_dictionary",
        project_root=project_root,
        phase="01_project_setup_and_governance",
    )

    try:
        summary = compile_data_dictionary(input_path, output_path, project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary, output=str(output_path))
    logger.write(status="ok")
    print(
        f"[OK]   Documented {summary['tables_documented']} table(s), {summary['total_fields']} field(s), in {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
