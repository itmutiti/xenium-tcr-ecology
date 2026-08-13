"""Append-only inventory CSV writer.

Header written once, paths relativised to the project root so
inventories are portable across machines, unknown fields rejected at
write time to keep each inventory's schema consistent.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from xenium_tcr_ecology.infra.exceptions import PipelineError


class InventoryWriter:
    """Appends rows to an audit-trail CSV/TSV, creating it with a header on
    first use. One row per tracked item (an excluded cell, a downloaded file,
    a probe audit result, ...) -- the exact schema is caller-defined via
    `fields`, matching whichever governance/reports table the phase script
    is producing (e.g. exclusion_log.tsv, patient_probe_audit.tsv)."""

    def __init__(self, path: Path, project_root: Path, fields: list[str], delimiter: str = "\t"):
        self.path = path
        self.project_root = project_root
        self.fields = fields
        self.delimiter = delimiter
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._is_new = not self.path.exists()

    def _relativize(self, value: Any) -> Any:
        if isinstance(value, Path):
            try:
                return str(value.relative_to(self.project_root))
            except ValueError:
                return str(value)
        return value

    def write_row(self, **fields: Any) -> None:
        unknown = set(fields) - set(self.fields)
        if unknown:
            raise PipelineError(f"Unknown inventory field(s): {sorted(unknown)}")
        row = {key: self._relativize(fields.get(key, "")) for key in self.fields}
        write_header = self._is_new
        with self.path.open("a", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=self.fields, delimiter=self.delimiter)
            if write_header:
                writer.writeheader()
                self._is_new = False
            writer.writerow(row)
