#!/usr/bin/env python3
"""
`13_clone_ecology_confirmatory_models/05_test_segmentation_robustness.py`

Segmentation robustness check for `04_quality_control/05_resegment_reference_subset.py`'s 3
representative resegmented sections -- see
src/xenium_tcr_ecology/clone_ecology/segmentation_robustness.py's module
docstring (including why a full resegmented-engagement recomputation is
explicitly out of this milestone's feasible scope).

Primary output: reports/clone_ecology/segmentation_robustness.pdf
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
from xenium_tcr_ecology.clone_ecology.segmentation_robustness import build_segmentation_robustness


def render_segmentation_robustness_report(
    comparison: pd.DataFrame, concordance: pd.DataFrame, output_path
) -> None:
    from xenium_tcr_ecology.viz.style import (
        COLORS,
        apply_publication_style,
        panel_label,
        panel_title,
    )

    apply_publication_style()
    # Side-by-side, not stacked: for 2 panels a single moderate-width
    # row keeps the aspect ratio closer to square than either a wide
    # 2-in-a-row or a tall 2x1 stack would.
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.4))

    y = np.arange(len(comparison))
    width = 0.32
    axes[0].barh(
        y - width / 2,
        comparison["mean_in_subset"],
        height=width,
        label="Resegmented subset (3 sections)",
        color=COLORS["sensitivity_analysis"],
    )
    axes[0].barh(
        y + width / 2,
        comparison["mean_outside_subset"],
        height=width,
        label="Rest of cohort",
        color=COLORS["primary_analysis"],
    )
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(comparison["metric"])
    axes[0].set_xlabel("Mean value (primary segmentation)")
    axes[0].tick_params(axis="both", which="major", pad=5)
    panel_title(axes[0], "Resegmented subset vs. cohort")
    # Legend placed below the axis, not inline: with only 3 categories
    # there is no empty quadrant inside the plot for an inline legend
    # without covering a bar.
    axes[0].legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1, handlelength=1.6, handletextpad=0.6
    )
    panel_label(axes[0], "A")

    axes[1].bar(
        concordance["section_id"],
        concordance["fraction_concordant_same_cell"],
        color=COLORS["sensitivity_analysis"],
        width=0.6,
    )
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Transcripts assigned to the same cell\n(primary vs. resegmented)")
    axes[1].tick_params(axis="both", which="major", pad=5)
    panel_title(axes[1], "Transcript-level concordance")
    # Default panel-letter x offset (-0.16) sits too close to this
    # panel's own long two-line y-axis label; shifted further left to
    # clear it.
    panel_label(axes[1], "B", x=-0.34)

    fig.tight_layout(w_pad=3.6, rect=(0, 0.1, 1, 0.92))
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
        logs_dir=project_root / "results" / "logs" / "13_clone_ecology_confirmatory_models",
        script_name="05_test_segmentation_robustness",
        project_root=project_root,
        phase="13_clone_ecology_confirmatory_models",
    )

    try:
        summary = build_segmentation_robustness(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    comparison = pd.read_parquet(summary["output_path"])
    concordance = pd.read_parquet(summary["concordance_output_path"])
    report_path = project_root / "reports" / "clone_ecology" / "segmentation_robustness.pdf"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    render_segmentation_robustness_report(comparison, concordance, report_path)
    summary["report_path"] = str(report_path)

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_resegmented_sections']} resegmented section(s), mean transcript concordance "
        f"{summary['mean_transcript_concordance']:.3f}. {summary['n_metrics_same_direction']}/{summary['n_metrics_total']} "
        f"engagement metrics and {summary['n_barrier_metrics_same_direction']}/{summary['n_barrier_metrics_total']} "
        f"barrier metrics (`14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`) same-sign in vs. outside the resegmented subset. Wrote {summary['output_path']}, "
        f"{summary['barrier_output_path']}, {summary['concordance_output_path']}, {summary['report_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
