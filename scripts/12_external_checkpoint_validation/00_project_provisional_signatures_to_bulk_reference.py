#!/usr/bin/env python3
"""
`12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`

Projects this project's T-cell-state marker methodology onto
GSE103322 (Puram et al. 2017), an independent, full-transcriptome
HNSCC scRNA-seq reference (data/external/GSE103322/README.md) -- see
src/xenium_tcr_ecology/external_checkpoint/bulk_reference.py's module docstring
.

Primary output: reports/external_checkpoint/bulk_projection.pdf
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.external_checkpoint.bulk_reference import build_bulk_projection


def render_bulk_projection_report(comparison: pd.DataFrame, output_path) -> None:
    from xenium_tcr_ecology.viz.style import COLORS, apply_publication_style, panel_title

    apply_publication_style()
    comparison_sorted = comparison.sort_values("project_fraction", ascending=True)
    y = np.arange(len(comparison_sorted))

    # Narrower, taller canvas: a word processor auto-scales a pasted
    # image to fit the page's text width (~6.5in), shrinking every
    # dimension including font size by that ratio. A narrower canvas
    # shrinks less, so the same in-image point size reads larger once
    # pasted.
    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    ax.hlines(
        y,
        comparison_sorted["project_fraction"],
        comparison_sorted["reference_fraction"],
        color=COLORS["reference_line"],
        linewidth=1.3,
        zorder=1,
    )
    ax.scatter(
        comparison_sorted["project_fraction"],
        y,
        color=COLORS["primary_analysis"],
        s=60,
        zorder=2,
        edgecolors="none",
        label="This study (Xenium)",
    )
    ax.scatter(
        comparison_sorted["reference_fraction"],
        y,
        color=COLORS["sensitivity_analysis"],
        s=60,
        zorder=2,
        marker="D",
        edgecolors="none",
        label="GSE103322 (Puram et al. 2017, external)",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(comparison_sorted["state"])
    ax.set_xlabel("Fraction of T cells")
    ax.tick_params(axis="both", which="major", pad=5)
    panel_title(ax, "T-cell state proportions vs. external reference")
    ax.legend(loc="lower right", handlelength=1.6, handletextpad=0.6)
    fig.tight_layout()
    fig.savefig(output_path)
    fig.savefig(str(output_path).replace(".pdf", ".png"), dpi=600)
    plt.close(fig)


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "12_external_checkpoint_validation",
        script_name="00_project_provisional_signatures_to_bulk_reference",
        project_root=project_root,
        phase="12_external_checkpoint_validation",
    )

    try:
        summary = build_bulk_projection(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    comparison = pd.read_parquet(summary["comparison_path"])
    report_path = project_root / "reports" / "external_checkpoint" / "bulk_projection.pdf"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    render_bulk_projection_report(comparison, report_path)
    summary["report_path"] = str(report_path)

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_reference_t_cells']} reference T cell(s) (GSE103322), "
        f"{summary['n_project_t_cells']} project T cell(s). Reference state counts: {summary['reference_state_counts']}. "
        f"Wrote {summary['output_path']}, {summary['comparison_path']}, {summary['report_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
