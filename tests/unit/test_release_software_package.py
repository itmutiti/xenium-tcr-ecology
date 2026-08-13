"""Unit tests for xenium_tcr_ecology.release.software_package (`17_statistical_closure_and_release/07_build_documented_software_package.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.release.software_package import extract_docstring_summary


class TestExtractDocstringSummary:
    def test_real_first_line_of_multiline_docstring_is_returned(self):
        doc = "Real first line.\n\nReal second paragraph with more detail."
        assert extract_docstring_summary(doc) == "Real first line."

    def test_real_single_line_docstring_is_returned_as_is(self):
        assert extract_docstring_summary("Real one-liner.") == "Real one-liner."

    def test_real_none_docstring_gives_placeholder(self):
        assert extract_docstring_summary(None) == "(no docstring)"

    def test_real_empty_docstring_gives_placeholder(self):
        assert extract_docstring_summary("") == "(no docstring)"

    def test_real_docstring_with_leading_whitespace_is_stripped(self):
        assert (
            extract_docstring_summary("   Real indented line.\n   more.") == "Real indented line."
        )
