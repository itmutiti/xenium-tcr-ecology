"""Unit tests for xenium_tcr_ecology.release.manuscript (`17_statistical_closure_and_release/05_render_reproducible_manuscript.qmd`)."""

from __future__ import annotations

from xenium_tcr_ecology.release.manuscript import strip_qmd_frontmatter_and_comments


class TestStripQmdFrontmatterAndComments:
    def test_real_yaml_frontmatter_is_removed(self):
        qmd = "---\ntitle: X\n---\n\nReal body text.\n"
        result = strip_qmd_frontmatter_and_comments(qmd)
        assert "title: X" not in result
        assert "Real body text." in result

    def test_real_html_comment_block_is_removed(self):
        qmd = "---\ntitle: X\n---\n\n<!--\nreal internal note\n-->\n\nReal body text.\n"
        result = strip_qmd_frontmatter_and_comments(qmd)
        assert "real internal note" not in result
        assert "Real body text." in result
