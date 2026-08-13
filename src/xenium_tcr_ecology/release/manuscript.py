"""Renders the reproducible-manuscript source
(`scripts/17_statistical_closure_and_release/05_render_reproducible_manuscript.qmd`)
to PDF via `markdown` + `xhtml2pdf`: pure Python, no system-level C
library dependencies (unlike `weasyprint`)."""

from __future__ import annotations

import re
from pathlib import Path

import markdown
from xhtml2pdf import pisa

from xenium_tcr_ecology.infra.exceptions import PipelineError

QMD_SOURCE_PATH = "scripts/17_statistical_closure_and_release/05_render_reproducible_manuscript.qmd"


def strip_qmd_frontmatter_and_comments(qmd_text: str) -> str:
    """Pure, testable: removes the YAML frontmatter (between the first
    two `---` lines) and any HTML comment blocks (`<!-- ... -->`),
    leaving only the Markdown body."""
    without_frontmatter = re.sub(r"^---\n.*?\n---\n", "", qmd_text, count=1, flags=re.DOTALL)
    without_comments = re.sub(r"<!--.*?-->", "", without_frontmatter, flags=re.DOTALL)
    return without_comments.strip()


def build_manuscript_pdf(project_root: Path) -> dict:
    qmd_path = project_root / QMD_SOURCE_PATH
    if not qmd_path.is_file():
        raise PipelineError(f"'{qmd_path}' not found.")

    qmd_text = qmd_path.read_text()
    body = strip_qmd_frontmatter_and_comments(qmd_text)
    if not body:
        raise PipelineError(
            f"'{qmd_path}' has no Markdown body after stripping frontmatter/comments -- not yet written."
        )

    html_body = markdown.markdown(body, extensions=["extra", "toc"])
    html_document = f"<html><head><meta charset='utf-8'></head><body>{html_body}</body></html>"

    manuscript_dir = project_root / "manuscript"
    manuscript_dir.mkdir(parents=True, exist_ok=True)
    output_path = manuscript_dir / "manuscript.pdf"

    with open(output_path, "wb") as f:
        result = pisa.CreatePDF(html_document, dest=f)
    if result.err:
        raise PipelineError(
            f"xhtml2pdf reported {result.err} error(s) rendering '{qmd_path}' to PDF."
        )

    return {
        "source_path": str(qmd_path),
        "output_path": str(output_path),
        "n_body_characters": len(body),
    }
