#!/usr/bin/env python3
"""
Renders a Markdown file to PDF -- the general-purpose tool behind
`docs/execution_manual/EXECUTION_MANUAL.pdf`'s build. The Markdown
source is always authoritative; a generated PDF is a convenience copy
for offline/print use, not an independent source -- if the two ever
disagree, re-run this script rather than hand-editing the PDF.

Uses the same tooling as the manuscript renderer
(`src/xenium_tcr_ecology/release/manuscript.py`): `markdown` +
`xhtml2pdf`, pure Python with no system C library dependencies.

Usage:
    python3 tools/build_markdown_pdf.py <input.md> [output.pdf]
    (output.pdf defaults to <input>.pdf next to the source)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import markdown
from xhtml2pdf import pisa

PDF_CSS = """
<style>
body { font-family: Helvetica, Arial, sans-serif; font-size: 10pt; line-height: 1.4; }
h1 { font-size: 18pt; margin-top: 0; }
h2 { font-size: 14pt; margin-top: 18pt; border-bottom: 1px solid #999; }
h3 { font-size: 12pt; margin-top: 14pt; }
code { font-family: Courier, monospace; background-color: #f0f0f0; padding: 1px 3px; }
pre { font-family: Courier, monospace; background-color: #f0f0f0; padding: 6px; font-size: 8pt; }
table { border-collapse: collapse; width: 100%; margin: 8pt 0; table-layout: fixed; }
th, td { border: 1px solid #999; padding: 4px 6px; font-size: 8.5pt; text-align: left;
         word-wrap: break-word; overflow-wrap: break-word; word-break: break-all; }
</style>
"""

# xhtml2pdf/reportlab does not enforce a fixed table-cell width against a
# long run of characters with no whitespace at all (e.g.
# "spatial_graph_construction_and_calibration", a single 43-character
# token) -- confirmed by rendering several candidate fixes before
# picking this one: `word-break: break-all` is not honoured, a zero-width
# space (U+200B) renders as a visible tofu box in this renderer's font
# handling rather than an invisible break, and even an explicit HTML
# <colgroup> width or a <wbr> break *opportunity* is ignored -- the cell
# simply expands into its neighbour regardless. Only a forced, literal
# `<br>` at chosen points reliably constrains it. So this inserts real
# `<br>` tags into long unbroken tokens, breaking at the nearest natural
# boundary (_, -, /, .) to a target line width -- this does change the
# rendered (and copy-pasted) text by adding line breaks, unlike a
# zero-width space, but nothing else tested actually works with this
# rendering engine.
#
# Two different thresholds, confirmed necessary by rendering both ways:
# ordinary paragraph/list flow *can* wrap a whole long token onto its own
# line without breaking it internally (normal reflow), and only overflows
# the page margin when a single token is long enough to exceed the full
# page width by itself -- applying the narrow table-cell threshold there
# too over-broke short, perfectly-fine-inline tokens like
# `environment/conda/environment.lock` into several stacked lines for no
# reason. A table *column* is much narrower than the full page, so it
# needs the aggressive threshold to avoid overflowing into the next
# column.
_LONG_TOKEN_RE_TABLE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_\-/.]{14,}")
_LONG_TOKEN_RE_BODY = re.compile(r"[A-Za-z0-9][A-Za-z0-9_\-/.]{49,}")
_TARGET_LINE_WIDTH_TABLE = 16
# Body text has the full page width available, not a narrow table column --
# confirmed by rendering: reusing the table's 16-char width here fragmented
# long paths into many tiny stacked lines even though they'd have fit on one
# or two normal-width lines. ~65 chars comfortably fits this document's body
# column at its font size without risking overflow for the longest tokens
# actually present.
_TARGET_LINE_WIDTH_BODY = 65


def _break_long_token(token: str, target_width: int) -> str:
    """Pure, testable: inserts <br> into a single long identifier-like
    token at natural boundaries, greedily packing each line up to
    ~target_width characters. A hyphen stays attached to the end of
    the line before it (conventional hyphenation, e.g. "transcript-" /
    "integrity"); underscore/slash/period stay attached to the start of
    the line after (path/identifier convention, e.g. "spatial_" /
    "graph"), so a break point never falls mid-word either way."""
    parts = re.split(r"([_\-/.])", token)
    pieces: list[str] = []
    buf = parts[0]
    for i in range(1, len(parts), 2):
        delim = parts[i]
        nxt = parts[i + 1] if i + 1 < len(parts) else ""
        if delim == "-":
            pieces.append(buf + delim)
            buf = nxt
        else:
            pieces.append(buf)
            buf = delim + nxt
    if buf:
        pieces.append(buf)
    pieces = [p for p in pieces if p]

    lines: list[str] = []
    current = ""
    for piece in pieces:
        if current and len(current) + len(piece) > target_width:
            lines.append(current)
            current = piece
        else:
            current += piece
    if current:
        lines.append(current)
    return "<br/>".join(lines)


def _add_soft_breaks(html_body: str) -> str:
    """Pure, testable: forces line breaks into long unbroken identifiers
    (file paths, `code_spans`, snake_case names) wherever they could
    overflow -- the aggressive threshold inside <td>/<th> (a table column
    is narrow), a much higher one everywhere else (ordinary paragraph/list
    flow already wraps a whole long-but-not-page-width token onto its own
    line; only a token long enough to exceed the full page width by
    itself needs breaking there). Skips <pre> blocks (multi-line code
    examples), where inserting a break would corrupt the literal command
    shown, and where the raw <pre> already wraps or scrolls acceptably."""

    def cell_repl(match: re.Match) -> str:
        open_tag, content, close_tag = match.group(1), match.group(2), match.group(3)
        content = _LONG_TOKEN_RE_TABLE.sub(
            lambda m: _break_long_token(m.group(0), _TARGET_LINE_WIDTH_TABLE), content
        )
        return open_tag + content + close_tag

    parts = re.split(r"(<pre>.*?</pre>)", html_body, flags=re.DOTALL)
    for i, part in enumerate(parts):
        if part.startswith("<pre>"):
            continue
        part = re.sub(r"(<t[dh][^>]*>)(.*?)(</t[dh]>)", cell_repl, part, flags=re.DOTALL)
        part = _LONG_TOKEN_RE_BODY.sub(
            lambda m: _break_long_token(m.group(0), _TARGET_LINE_WIDTH_BODY), part
        )
        parts[i] = part
    return "".join(parts)


def build_markdown_pdf(source_path: Path, output_path: Path) -> dict:
    if not source_path.is_file():
        raise FileNotFoundError(f"'{source_path}' not found.")

    md_text = source_path.read_text()
    html_body = markdown.markdown(md_text, extensions=["extra", "toc"])
    html_body = _add_soft_breaks(html_body)
    html_document = (
        f"<html><head><meta charset='utf-8'>{PDF_CSS}</head><body>{html_body}</body></html>"
    )

    with open(output_path, "wb") as f:
        result = pisa.CreatePDF(html_document, dest=f)
    if result.err:
        raise RuntimeError(f"xhtml2pdf reported {result.err} error(s) rendering '{source_path}'.")

    return {
        "source_path": str(source_path),
        "output_path": str(output_path),
        "n_source_characters": len(md_text),
    }


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python3 tools/build_markdown_pdf.py <input.md> [output.pdf]", file=sys.stderr)
        return 1

    source_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else source_path.with_suffix(".pdf")

    try:
        summary = build_markdown_pdf(source_path, output_path)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    print(
        f"[OK]   Rendered {summary['source_path']} -> {summary['output_path']} ({summary['n_source_characters']} source characters)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
