"""Unit tests for xenium_tcr_ecology.release.main_figures (`17_statistical_closure_and_release/02_generate_all_main_figures.py`)."""

from __future__ import annotations

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.release.main_figures import MAIN_FIGURE_MANIFEST, build_main_figures


def _real_minimal_pdf_bytes() -> bytes:
    """A real, valid, single-page PDF (not placeholder text) -- needed
    because `build_main_figures` now also converts each PDF to PNG/SVG
    via `pdftocairo`, which requires real PDF structure, not arbitrary
    bytes with a `%PDF` prefix."""
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3])
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    plt.close(fig)
    return buf.getvalue()


def _touch_source_pdfs(project_root):
    pdf_bytes = _real_minimal_pdf_bytes()
    for entry in MAIN_FIGURE_MANIFEST:
        path = project_root / entry["source_path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(pdf_bytes)


class TestBuildMainFigures:
    def test_real_all_figures_copied_with_sequential_numbering(self, tmp_path):
        _touch_source_pdfs(tmp_path)
        summary = build_main_figures(tmp_path)
        assert summary["n_figures"] == len(MAIN_FIGURE_MANIFEST)
        figures_dir = tmp_path / "figures" / "main"
        first = MAIN_FIGURE_MANIFEST[0]
        assert (
            figures_dir / f"Figure_{first['figure_number']}_{first['descriptive_name']}.pdf"
        ).exists()

    def test_real_png_and_svg_also_exported(self, tmp_path):
        _touch_source_pdfs(tmp_path)
        build_main_figures(tmp_path)
        figures_dir = tmp_path / "figures" / "main"
        first = MAIN_FIGURE_MANIFEST[0]
        stem = f"Figure_{first['figure_number']}_{first['descriptive_name']}"
        assert (figures_dir / f"{stem}.png").exists()
        assert (figures_dir / f"{stem}.svg").exists()

    def test_real_missing_source_pdf_raises(self, tmp_path):
        with pytest.raises(PipelineError):
            build_main_figures(tmp_path)

    def test_real_stale_figure_from_a_previous_manifest_is_removed(self, tmp_path):
        _touch_source_pdfs(tmp_path)
        figures_dir = tmp_path / "figures" / "main"
        figures_dir.mkdir(parents=True, exist_ok=True)
        stale = figures_dir / "Figure_1_Some_Superseded_Descriptive_Name.pdf"
        stale.write_bytes(_real_minimal_pdf_bytes())
        build_main_figures(tmp_path)
        assert not stale.exists()
