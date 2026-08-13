"""Assembles the complete set of supplementary figures -- every
diagnostic/sensitivity/robustness PDF this project generated across all
17 computational workflow stages that is not already used as a main
figure (`17_statistical_closure_and_release/02_generate_all_main_figures.py`), plus the complete clone atlas
(`13_clone_ecology_confirmatory_models/06_generate_clone_atlas.py`, all 261
clone cards -- the "non-selected clone maps" this milestone's scaffold
names) -- `17_statistical_closure_and_release/03_generate_all_supplementary_figures.py`.

**Deliberate completeness, not a hand-picked subset:** every report PDF
this project has ever generated (28 PDFs at time of writing, across QC,
preprocessing, tumour, TCR, graphs, niches, clone ecology, external
checkpoint, interactions, HPV and validation) is included unless it is
one of the 6 already selected as a main figure (`17_statistical_closure_and_release/02_generate_all_main_figures.py`) --
matching this milestone's scaffold wording ("complete diagnostics,
sensitivity analyses..."), not a curated highlight reel. Each filename
uses the source report's already-descriptive stem (e.g.
`clone_size_sensitivity`, `null_model_calibration`) -- no phase/stage
numbers ever appeared in these names.

**Multi-format export, not PDF-only:** every figure PDF is also
rendered to a high-resolution PNG (300 DPI) and an SVG alongside it
(`xenium_tcr_ecology.infra.figure_export`) -- see `main_figures.py`'s
module docstring for the full rationale. The clone atlas
(`SupplementaryClone_Atlas.html`) is not a figure PDF and is not
converted.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.figure_export import export_pdf_to_png_and_svg
from xenium_tcr_ecology.release.main_figures import MAIN_FIGURE_MANIFEST

CLONE_ATLAS_PATH = "reports/clone_ecology/clone_atlas.html"


def find_supplementary_pdfs(project_root: Path, main_figure_sources: set[str]) -> list[str]:
    """Pure, testable: sorted, deterministic list of every `.pdf` under
    `reports/` (relative to `project_root`), excluding
    `main_figure_sources`."""
    reports_dir = project_root / "reports"
    if not reports_dir.is_dir():
        return []
    all_pdfs = sorted(str(p.relative_to(project_root)) for p in reports_dir.rglob("*.pdf"))
    return [p for p in all_pdfs if p not in main_figure_sources]


def build_supplementary_figures(project_root: Path) -> dict:
    figures_dir = project_root / "figures" / "supplementary"
    main_figure_sources = {entry["source_path"] for entry in MAIN_FIGURE_MANIFEST}

    supplementary_pdfs = find_supplementary_pdfs(project_root, main_figure_sources)
    if not supplementary_pdfs:
        raise PipelineError(f"No supplementary PDFs found under '{project_root / 'reports'}'.")

    clone_atlas_path = project_root / CLONE_ATLAS_PATH
    if not clone_atlas_path.is_file():
        raise PipelineError(
            f"'{clone_atlas_path}' not found. Run `13_clone_ecology_confirmatory_models/06_generate_clone_atlas.py` first."
        )

    figures_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for i, rel_path in enumerate(supplementary_pdfs, start=1):
        source = project_root / rel_path
        stem = f"SupplementaryFigure_{i}_{Path(rel_path).stem}"
        dest = figures_dir / f"{stem}.pdf"
        dest.write_bytes(source.read_bytes())
        export_pdf_to_png_and_svg(dest, figures_dir, stem)
        manifest_rows.append(
            {"supplementary_figure_number": i, "filename": dest.name, "source_path": rel_path}
        )

    atlas_dest = figures_dir / "SupplementaryClone_Atlas.html"
    atlas_dest.write_bytes(clone_atlas_path.read_bytes())
    manifest_rows.append(
        {
            "supplementary_figure_number": len(supplementary_pdfs) + 1,
            "filename": atlas_dest.name,
            "source_path": CLONE_ATLAS_PATH,
        }
    )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = figures_dir / "MANIFEST.tsv"
    manifest_df.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_supplementary_figures": len(manifest_rows),
        "figures_dir": str(figures_dir),
        "manifest_path": str(manifest_path),
    }
