"""Unit tests for xenium_tcr_ecology.release.supplementary_figures (`17_statistical_closure_and_release/03_generate_all_supplementary_figures.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.release.supplementary_figures import find_supplementary_pdfs


class TestFindSupplementaryPdfs:
    def test_real_excludes_main_figure_sources(self, tmp_path):
        reports_dir = tmp_path / "reports"
        (reports_dir / "a").mkdir(parents=True)
        (reports_dir / "a" / "included.pdf").write_bytes(b"%PDF real")
        (reports_dir / "a" / "excluded.pdf").write_bytes(b"%PDF real")
        main_figure_sources = {"reports/a/excluded.pdf"}
        result = find_supplementary_pdfs(tmp_path, main_figure_sources)
        assert result == ["reports/a/included.pdf"]

    def test_real_returns_sorted_deterministic_order(self, tmp_path):
        reports_dir = tmp_path / "reports"
        (reports_dir / "z").mkdir(parents=True)
        (reports_dir / "a").mkdir(parents=True)
        (reports_dir / "z" / "z.pdf").write_bytes(b"%PDF real")
        (reports_dir / "a" / "a.pdf").write_bytes(b"%PDF real")
        result = find_supplementary_pdfs(tmp_path, set())
        assert result == sorted(result)
