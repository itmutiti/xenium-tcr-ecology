"""Multi-format figure export: every assembled manuscript/supplementary
figure is a vector PDF (publication-quality, kept as the primary
format) plus a high-resolution PNG (600 DPI, meeting typical journal
print-reproduction minimums, for easy inclusion in slide decks and
Word/PowerPoint, which do not embed PDF well), an SVG (vector, for
direct editing in Illustrator/Inkscape or embedding in HTML), and a
600 DPI TIFF (LZW-compressed, lossless -- some submission portals'
separate-figure-upload paths accept only JPEG/TIFF/EPS, not PDF or
PNG; TIFF is the one of those three that stays lossless). No JPEG --
lossy compression is inappropriate for statistical plots with sharp
text and thin lines.

Uses `pdftocairo` (poppler-utils, already a system dependency of this
environment) to convert the already-generated, already-reviewed
PDF directly, rather than re-plotting from each analysis script in a
second format: this keeps the ~30 individual plotting scripts focused on
producing one correct PDF, and centralises multi-format export in one
place, applied uniformly to every assembled figure.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from xenium_tcr_ecology.infra.exceptions import PipelineError

PNG_DPI = 600
TIFF_DPI = 600


def check_pdftocairo_available() -> None:
    if shutil.which("pdftocairo") is None:
        raise PipelineError(
            "'pdftocairo' (poppler-utils) not found on PATH -- required for PNG/SVG/TIFF figure "
            "export. Install poppler-utils (e.g. `apt-get install poppler-utils` / "
            "`conda install poppler`)."
        )


def export_pdf_to_png_svg_and_tiff(pdf_path: Path, output_dir: Path, stem: str) -> dict:
    """Renders `pdf_path` (a single-page figure PDF) to a PNG
    (`PNG_DPI` dots per inch), an SVG, and an LZW-compressed TIFF
    (`TIFF_DPI` dots per inch), all named `<stem>` in `output_dir`.
    Returns their paths."""
    check_pdftocairo_available()

    png_path = output_dir / f"{stem}.png"
    svg_path = output_dir / f"{stem}.svg"
    tiff_path = output_dir / f"{stem}.tif"

    subprocess.run(
        [
            "pdftocairo",
            "-png",
            "-r",
            str(PNG_DPI),
            "-singlefile",
            str(pdf_path),
            str(output_dir / stem),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["pdftocairo", "-svg", str(pdf_path), str(svg_path)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            "pdftocairo",
            "-tiff",
            "-tiffcompression",
            "lzw",
            "-r",
            str(TIFF_DPI),
            "-singlefile",
            str(pdf_path),
            str(output_dir / stem),
        ],
        check=True,
        capture_output=True,
    )

    if not png_path.is_file() or not svg_path.is_file() or not tiff_path.is_file():
        raise PipelineError(
            f"pdftocairo did not produce expected PNG/SVG/TIFF output for '{pdf_path}' "
            f"(stem '{stem}')."
        )

    return {"png_path": str(png_path), "svg_path": str(svg_path), "tiff_path": str(tiff_path)}
