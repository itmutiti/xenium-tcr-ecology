#!/usr/bin/env python3
"""
`16_external_validation_and_generalisation/04_validate_ecosystem_signatures_in_bulk.py`

Projects Niche and Ecosystem Discovery ecosystem-derived gene signatures onto the
TCGA-HNSC bulk RNA-seq cohort and cautiously tests correlation with
PTPRC (CD45), an independent immune-infiltration proxy. See
src/xenium_tcr_ecology/validation/bulk_projection.py's module docstring
.

Primary output: reports/validation/bulk_projection.pdf
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.validation.bulk_projection import build_bulk_projection


def render_bulk_projection_report(result: pd.DataFrame, output_path) -> None:
    from xenium_tcr_ecology.viz.style import (
        COLORS,
        FS_LEGEND,
        OKABE_ITO,
        apply_publication_style,
        panel_title,
    )

    apply_publication_style()
    # Narrower, taller canvas: a word processor auto-scales a pasted
    # image to fit the page's text width (~6.5in), shrinking every
    # dimension including font size by that ratio. A narrower canvas
    # shrinks less, so the same in-image point size reads larger once
    # pasted.
    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    colors = [
        (
            OKABE_ITO["bluish_green"]
            if m
            else (OKABE_ITO["vermillion"] if m is False else COLORS["not_significant"])
        )
        for m in result["matches_expected_direction"]
    ]
    y = range(len(result))
    ax.hlines(y, 0, result["rho_vs_immune_proxy"], color=colors, linewidth=6, zorder=1)
    ax.set_yticks(list(y))
    ax.set_yticklabels(result["ecosystem_label"])
    ax.axvline(0, linestyle="--", color=COLORS["reference_line"], linewidth=1.1, zorder=0)
    ax.set_xlabel("Spearman rho vs. PTPRC (CD45) expression")
    ax.tick_params(axis="both", which="major", pad=5)
    panel_title(ax, "Ecosystem signature vs. immune-infiltration proxy (TCGA-HNSC)")
    handles = [
        plt.Line2D(
            [0],
            [0],
            color=OKABE_ITO["bluish_green"],
            linewidth=6,
            label="Matches predicted direction",
        ),
        plt.Line2D([0], [0], color=OKABE_ITO["vermillion"], linewidth=6, label="Does not match"),
        plt.Line2D(
            [0],
            [0],
            color=COLORS["not_significant"],
            linewidth=6,
            label="No pre-specified prediction",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        fontsize=FS_LEGEND - 1,
        handlelength=1.8,
        handletextpad=0.6,
    )
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
        logs_dir=project_root / "results" / "logs" / "16_external_validation_and_generalisation",
        script_name="04_validate_ecosystem_signatures_in_bulk",
        project_root=project_root,
        phase="16_external_validation_and_generalisation",
    )

    try:
        summary = build_bulk_projection(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    result = pd.read_parquet(summary["output_path"])
    report_path = project_root / "reports" / "validation" / "bulk_projection.pdf"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    render_bulk_projection_report(result, report_path)
    summary["report_path"] = str(report_path)

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_tcga_samples']} TCGA-HNSC sample(s), {summary['n_ecosystems_tested']} ecosystem signature(s) tested. "
        f"{summary['n_matching_expected_direction']}/{summary['n_directional_predictions_made']} match the expected direction. "
        f"Wrote {summary['output_path']}, {report_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
