#!/usr/bin/env python3
"""
`10_niche_and_ecosystem_discovery/07_leave_one_patient_out_niche_stability.py`

For each of `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s cohort patients, refits archetype
clustering with that patient withheld entirely and checks (1) whether
the refit centroids still resemble the full-dataset archetypes, and (2)
whether the withheld patient's own cells are still correctly identified
from the remaining patients' structure -- see
src/xenium_tcr_ecology/niches/lopo_stability.py's module docstring.

Primary output: reports/niches/LOPO_stability.pdf
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
from xenium_tcr_ecology.niches.lopo_stability import build_lopo_stability


def render_lopo_stability_report(result: pd.DataFrame, output_path) -> None:
    from xenium_tcr_ecology.viz.style import (
        COLORS,
        apply_publication_style,
        panel_label,
        panel_title,
    )

    apply_publication_style()
    per_patient = result.drop_duplicates("patient_id").sort_values("identifiability_accuracy")

    # Side-by-side, not stacked: for 2 panels a single moderate-width
    # row keeps the aspect ratio closer to square than either a wide
    # 2-in-a-row or a tall 2x1 stack would.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.8))

    ax1.barh(
        per_patient["patient_id"],
        per_patient["identifiability_accuracy"],
        color=COLORS["primary_analysis"],
        height=0.65,
    )
    ax1.axvline(
        1.0 / result["archetype"].nunique(),
        linestyle="--",
        color=COLORS["reference_line"],
        linewidth=1.1,
        label="Chance level (1/K)",
    )
    ax1.set_xlabel("Identifiability accuracy\n(withheld patient's own cells)")
    ax1.set_xlim(0, 1.0)
    ax1.tick_params(axis="both", which="major", pad=5)
    # Legend placed below the axis, not inline: with 11 patient bars
    # spanning most of the plot area, there is no reliably empty
    # quadrant for an inline legend without covering a bar.
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), handlelength=1.8, handletextpad=0.6)
    panel_title(ax1, "Withheld-patient identifiability")
    panel_label(ax1, "A")

    archetype_order = sorted(result["archetype"].unique())
    data_by_archetype = [
        result.loc[result["archetype"] == a, "centroid_cosine_similarity"] for a in archetype_order
    ]
    bp = ax2.boxplot(
        data_by_archetype,
        tick_labels=[str(a) for a in archetype_order],
        patch_artist=True,
        widths=0.5,
    )
    for box in bp["boxes"]:
        box.set(facecolor="white", edgecolor=COLORS["not_significant"], linewidth=1.4)
    for median in bp["medians"]:
        median.set(color=COLORS["not_significant"], linewidth=1.7)
    for i, values in enumerate(data_by_archetype, start=1):
        ax2.scatter(
            [i] * len(values),
            values,
            alpha=0.7,
            color=COLORS["sensitivity_analysis"],
            s=32,
            zorder=3,
            edgecolors="none",
        )
    ax2.set_xlabel("Archetype")
    ax2.set_ylabel("Centroid cosine similarity\n(LOPO refit vs. full-dataset centroid)")
    ax2.set_ylim(0, 1.05)
    ax2.tick_params(axis="both", which="major", pad=5)
    panel_title(ax2, "Centroid stability across LOPO refits")
    # Default panel-letter x offset (-0.16) sits too close to this
    # panel's own long two-line y-axis label; shifted further left to
    # clear it.
    panel_label(ax2, "B", x=-0.30)

    fig.tight_layout(w_pad=3.6, rect=(0, 0.1, 1, 1))
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
        logs_dir=project_root / "results" / "logs" / "10_niche_and_ecosystem_discovery",
        script_name="07_leave_one_patient_out_niche_stability",
        project_root=project_root,
        phase="10_niche_and_ecosystem_discovery",
    )

    try:
        summary = build_lopo_stability(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    result = pd.read_parquet(summary["output_path"])
    report_path = project_root / "reports" / "niches" / "LOPO_stability.pdf"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    render_lopo_stability_report(result, report_path)
    summary["report_path"] = str(report_path)

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_patients']} patient(s), {summary['n_archetypes']} archetype(s). "
        f"Mean identifiability accuracy {summary['mean_identifiability_accuracy']:.3f} "
        f"(min {summary['min_identifiability_accuracy']:.3f}). Mean centroid cosine similarity "
        f"{summary['mean_centroid_cosine_similarity']:.3f} (min {summary['min_centroid_cosine_similarity']:.3f}). "
        f"Wrote {summary['output_path']}, {summary['report_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
