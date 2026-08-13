#!/usr/bin/env python3
"""
`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`

Confirms each planned null model (constrained permutation, degree-
preserving, graph-preserving -- the registered `null_model_calibration_suite`
analysis, governance/analysis_registry.tsv) achieves nominal Type I error
and adequate power at n=10-patient scale on `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s synthetic
ground-truth data, before being trusted on real clones. See
src/xenium_tcr_ecology/graphs/null_models.py's module docstring for the
three null model definitions.

Primary output: reports/graphs/null_model_calibration.pdf
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.graphs.null_models import NOMINAL_ALPHA, build_null_model_calibration

MODEL_LABELS = {
    "pvalue_constrained_permutation": "Constrained permutation",
    "pvalue_degree_preserving": "Degree-preserving",
    "pvalue_graph_preserving": "Graph-preserving",
}


MODEL_COLOR_KEYS = {
    "pvalue_constrained_permutation": "null_constrained_permutation",
    "pvalue_degree_preserving": "null_degree_preserving",
    "pvalue_graph_preserving": "null_graph_preserving",
}
MODEL_MARKERS = {
    "pvalue_constrained_permutation": "o",
    "pvalue_degree_preserving": "s",
    "pvalue_graph_preserving": "^",
}


def render_calibration_plot(calibration_summary: dict, output_path) -> None:
    from xenium_tcr_ecology.viz.style import (
        COLORS,
        apply_publication_style,
        panel_title,
        significance_reference_line,
    )

    apply_publication_style()
    # Narrower, taller canvas: a word processor auto-scales a pasted
    # image to fit the page's text width (~6.5in), shrinking every
    # dimension including font size by that ratio. A narrower canvas
    # shrinks less, so the same in-image point size reads larger once
    # pasted.
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    for col, by_effect in calibration_summary.items():
        sorted_keys = sorted(by_effect, key=float)
        effect_sizes = [float(k) for k in sorted_keys]
        rates = [by_effect[k]["rejection_rate"] for k in sorted_keys]
        ax.plot(
            effect_sizes,
            rates,
            marker=MODEL_MARKERS.get(col, "o"),
            markersize=6.5,
            markeredgewidth=0,
            linewidth=1.9,
            color=COLORS.get(MODEL_COLOR_KEYS.get(col), None),
            label=MODEL_LABELS.get(col, col),
        )
    significance_reference_line(ax, y=NOMINAL_ALPHA)
    ax.set_xlabel("Effect size (0 = true null)")
    ax.set_ylabel("Rejection rate at p < 0.05\n(10 synthetic patients)")
    ax.tick_params(axis="both", which="major", pad=5)
    panel_title(ax, "Null-model calibration: synthetic data")
    ax.legend(loc="lower right", handlelength=2.0, handletextpad=0.6)
    ax.set_ylim(-0.05, 1.05)
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
        logs_dir=project_root
        / "results"
        / "logs"
        / "09_spatial_graph_construction_and_calibration",
        script_name="08_run_calibration_suite_on_synthetic_data",
        project_root=project_root,
        phase="09_spatial_graph_construction_and_calibration",
    )

    try:
        summary = build_null_model_calibration(project_root)
        output_path = project_root / "reports" / "graphs" / "null_model_calibration.pdf"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        render_calibration_plot(summary["calibration_summary"], output_path)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_replicates']} replicate(s), {summary['n_permutations']} permutations each."
    )
    for model, by_effect in summary["calibration_summary"].items():
        print(f"       {model}:")
        for effect_size, stats in sorted(by_effect.items(), key=lambda kv: float(kv[0])):
            print(
                f"         effect={effect_size}: rejection_rate={stats['rejection_rate']} "
                f"({stats['n_rejected']}/{stats['n_total']}, 95% CI {stats['ci_95']})"
            )
    print(f"[OK]   Wrote {summary['output_path']}, {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
