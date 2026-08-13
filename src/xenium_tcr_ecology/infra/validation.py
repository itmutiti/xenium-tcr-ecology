"""Generic structured-input contract validation.

Fails loudly, before any downstream use, if a structured input is
missing required fields, has placeholder content left in, or falls short
of a minimum row count, rather than propagating an incomplete governance
record. Used by Literature and Novelty Audit's literature-audit scripts
and intended for reuse at every later stage boundary that consumes a
human-curated structured input.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from xenium_tcr_ecology.infra.exceptions import PipelineError

PLACEHOLDER_MARKERS = ("TODO", "TBD", "FIXME", "<TODO>")


def validate_records(
    records: list[dict[str, Any]],
    required_fields: list[str],
    source_path: Path,
    min_records: int = 1,
    no_placeholder_fields: list[str] | None = None,
) -> None:
    """Validate a list of dict records loaded from a structured input file.

    Raises PipelineError, listing every problem found (not just the first),
    if:
      - fewer than `min_records` records are present,
      - any record is missing a required field or has it empty,
      - any record's `no_placeholder_fields` still contain a placeholder
        marker (TODO/TBD/FIXME), indicating the human-curation step was not
        actually completed.
    """
    problems: list[str] = []

    if len(records) < min_records:
        problems.append(
            f"Expected at least {min_records} record(s) in '{source_path}', found {len(records)}."
        )

    no_placeholder_fields = no_placeholder_fields or []

    for i, record in enumerate(records):
        missing = [f for f in required_fields if not str(record.get(f, "")).strip()]
        if missing:
            problems.append(
                f"Record {i} in '{source_path}' is missing/empty required field(s): {missing}"
            )
        for field in no_placeholder_fields:
            value = str(record.get(field, ""))
            for marker in PLACEHOLDER_MARKERS:
                if marker in value:
                    problems.append(
                        f"Record {i} in '{source_path}', field '{field}' still contains "
                        f"placeholder marker '{marker}': {value!r}"
                    )

    if problems:
        raise PipelineError(
            f"Contract validation failed for '{source_path}' ({len(problems)} problem(s)):\n  "
            + "\n  ".join(problems)
        )
