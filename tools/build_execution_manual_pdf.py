#!/usr/bin/env python3
"""
Renders `docs/execution_manual/EXECUTION_MANUAL.md` (the single source of
truth) to `docs/execution_manual/EXECUTION_MANUAL.pdf` (a convenience copy
for offline/print use). The Markdown file is authoritative; if the two
ever disagree, re-run this script rather than hand-editing the PDF.

Thin wrapper around `tools/build_markdown_pdf.py` (the general-purpose
renderer) with this file's fixed source/output paths, kept as its own
script since it is the one
Markdown->PDF build referenced directly by name throughout `docs/
execution_manual/EXECUTION_MANUAL.md` and `README.md`.

Usage:
    python3 tools/build_execution_manual_pdf.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_markdown_pdf import build_markdown_pdf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "docs" / "execution_manual" / "EXECUTION_MANUAL.md"
OUTPUT_PATH = ROOT / "docs" / "execution_manual" / "EXECUTION_MANUAL.pdf"


def main() -> int:
    try:
        summary = build_markdown_pdf(SOURCE_PATH, OUTPUT_PATH)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(
        f"[OK]   Rendered {summary['source_path']} -> {summary['output_path']} ({summary['n_source_characters']} source characters)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
