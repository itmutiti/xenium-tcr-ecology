"""Redesigned, composite manuscript main figures
(`17_statistical_closure_and_release/11_build_redesigned_manuscript_figures.py`).

Builds multi-panel composite figures from already-computed, already-
validated data (no new statistical analysis is performed here),
replacing three of the original six main figures (`main_figures.py`) and
consolidating two others. See `docs/analysis_amendments.md` for why the
figure set was redesigned.

1. Framework generalisation (synthetic calibration + 2 independent
   tumour types) -- was Figure 1, still Figure 1, now 3 panels.
2. TCR probe validation against paired VDJ data -- new, previously had
   no figure.
3. Discrete vs. continuous clone structure -- unchanged, renumbered.
4. Variance partition with its feature-set sensitivity check shown
   directly -- was Figure 2, redesigned.
5. Barrier topology with its covariate-ablation robustness check and
   literature benchmark shown directly -- was Figure 4, redesigned.
6. HPV status, consolidated from two figures into one, with the
   multi-modal discordance QC finding included -- was Figures 5+6.

Every panel is built from an already-computed, already-governance-
documented data file.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from xenium_tcr_ecology.graphs.null_models import NOMINAL_ALPHA, clopper_pearson_ci
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.viz.style import (
    COLORS,
    FS_LEGEND,
    FS_PANEL_TITLE,
    apply_publication_style,
    panel_label,
    panel_title,
    significance_reference_line,
)

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
NULL_MODEL_COLS = list(MODEL_LABELS.keys())


def _calibration_summary(results: pd.DataFrame) -> dict:
    """Pure, testable: per-effect-size rejection rate + 95%
    Clopper-Pearson CI for each of the 3 null models -- the same
    computation `framework_generalisation.py` and `null_models.py`'s
    calibration-suite builder already use, factored out here so all
    three calibration panels (synthetic, breast, colorectal) use
    identical logic."""
    summary = {}
    for col in NULL_MODEL_COLS:
        by_effect = {}
        for effect_size, group in results.groupby("effect_size"):
            n_rejected = int((group[col] < NOMINAL_ALPHA).sum())
            n_total = len(group)
            ci = clopper_pearson_ci(n_rejected, n_total)
            by_effect[float(effect_size)] = {
                "rejection_rate": n_rejected / n_total,
                "ci_low": ci[0],
                "ci_high": ci[1],
            }
        summary[col] = by_effect
    return summary


MODEL_MARKERS = {
    "pvalue_constrained_permutation": "o",
    "pvalue_degree_preserving": "s",
    "pvalue_graph_preserving": "^",
}
# Small, deliberate horizontal offset per series so near-identical curves
# (the three null models converge almost everywhere) remain visually
# separable instead of one line fully occluding the other two.
MODEL_JITTER = {
    "pvalue_constrained_permutation": -0.006,
    "pvalue_degree_preserving": 0.0,
    "pvalue_graph_preserving": 0.006,
}


def _plot_calibration_panel(ax, summary: dict, title: str) -> None:
    for col in NULL_MODEL_COLS:
        by_effect = summary[col]
        effect_sizes = sorted(by_effect)
        x = np.array(effect_sizes) + MODEL_JITTER[col]
        rates = [by_effect[e]["rejection_rate"] for e in effect_sizes]
        ci_low = [by_effect[e]["ci_low"] for e in effect_sizes]
        ci_high = [by_effect[e]["ci_high"] for e in effect_sizes]
        color = COLORS[MODEL_COLOR_KEYS[col]]
        # Error bars drawn first (behind) so the heavier trend line and
        # markers sit visually on top, not competing with the CI for
        # attention.
        ax.errorbar(
            x,
            rates,
            yerr=[np.array(rates) - np.array(ci_low), np.array(ci_high) - np.array(rates)],
            fmt="none",
            ecolor=color,
            elinewidth=0.85,
            capsize=0,
            alpha=0.45,
            zorder=1,
        )
        ax.plot(
            x,
            rates,
            marker=MODEL_MARKERS[col],
            markersize=6.5,
            markeredgewidth=0,
            linewidth=1.9,
            color=color,
            label=MODEL_LABELS[col],
            zorder=2,
        )
    significance_reference_line(ax, y=NOMINAL_ALPHA)
    ax.set_xlabel("Effect size (0 = true null)")
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlim(-0.03, 1.03)
    ax.tick_params(axis="both", which="major", pad=5)
    panel_title(ax, title)


def build_figure_1_framework_generalisation(project_root: Path) -> dict:
    synthetic_path = project_root / "reports" / "graphs" / "null_model_calibration.parquet"
    breast_path = project_root / "data" / "derived" / "framework_generalisation_results.parquet"
    colorectal_path = (
        project_root
        / "data"
        / "derived"
        / "framework_generalisation_results_second_dataset.parquet"
    )
    output_path = (
        project_root
        / "reports"
        / "manuscript_figures"
        / "framework_generalisation_three_tumour_types.pdf"
    )

    for p in (synthetic_path, breast_path, colorectal_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    synthetic = _calibration_summary(pd.read_parquet(synthetic_path))
    breast = _calibration_summary(pd.read_parquet(breast_path))
    colorectal = _calibration_summary(pd.read_parquet(colorectal_path))

    apply_publication_style()
    # 2x2 grid, not a single row of 3 (too wide once shrunk to fit a
    # word-processor page) or a column of 3 (too tall/long once
    # rendered at a legible font size). A: top-left, B: top-right,
    # C: bottom-left, with the shared legend moved into the otherwise
    # empty bottom-right cell instead of a fourth row.
    # constrained_layout, not tight_layout(): tight_layout() ignores an
    # explicit GridSpec's own spacing and can silently re-compress rows
    # into each other (observed as panel-title/label overlap between
    # rows); constrained_layout is GridSpec-aware and resolves overlaps
    # itself.
    fig = plt.figure(figsize=(12.6, 10.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, hspace=0.15, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_legend = fig.add_subplot(gs[1, 1])
    ax_legend.axis("off")

    _plot_calibration_panel(ax_a, synthetic, "Synthetic calibration (n = 10)")
    _plot_calibration_panel(ax_b, breast, "Breast cancer (external)")
    _plot_calibration_panel(ax_c, colorectal, "Colorectal cancer (external)")
    for ax in (ax_a, ax_b, ax_c):
        ax.set_ylabel("Rejection rate at p < 0.05\n(95% Clopper-Pearson CI)")
    for letter, ax in zip("ABC", (ax_a, ax_b, ax_c)):
        panel_label(ax, letter)
    handles, labels = ax_a.get_legend_handles_labels()
    ax_legend.legend(
        handles,
        labels,
        loc="center",
        ncol=1,
        handlelength=2.4,
        columnspacing=2.4,
        handletextpad=0.7,
        frameon=False,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    fig.savefig(output_path.with_suffix(".png"), dpi=600)
    plt.close(fig)
    return {"output_path": str(output_path)}


def build_figure_variance_partition_with_sensitivity(project_root: Path) -> dict:
    original_path = project_root / "data" / "derived" / "variance_partition_results.parquet"
    sensitivity_path = (
        project_root
        / "data"
        / "derived"
        / "variance_partition_sensitivity_excluding_cycling.parquet"
    )
    output_path = (
        project_root / "reports" / "manuscript_figures" / "variance_partition_with_sensitivity.pdf"
    )

    for p in (original_path, sensitivity_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    original = pd.read_parquet(original_path).set_index("component")
    sensitivity = pd.read_parquet(sensitivity_path).set_index("component")

    # Ordered by the primary (11-feature) point estimate, largest first,
    # so the ranking the text reports (context > patient > identity) is
    # read directly off the panel from top to bottom.
    components = original.sort_values("proportion", ascending=False).index.tolist()
    component_labels = {
        "context": "Spatial context",
        "patient": "Patient",
        "identity": "Clonal identity",
    }

    apply_publication_style()
    # Narrower, taller canvas (was 9.6x3.9in): a word processor
    # auto-scales a pasted image to fit the page's text width (~6.5in),
    # shrinking every dimension including font size by that ratio. A
    # narrower canvas shrinks less, so the same in-image point size
    # reads larger once pasted. Here the point estimate AND its CI are
    # both primary content (a forest/point-range plot), unlike Figure
    # 1's calibration curves where the CI is secondary context -- so
    # error bars are kept fully opaque and made modestly thicker for
    # legibility, not lightened.
    fig, ax = plt.subplots(1, 1, figsize=(7.4, 5.6))
    y = np.arange(len(components))
    dodge = 0.14

    for offset, data, key, label in [
        (+dodge, original, "primary_analysis", "11 features (prespecified)"),
        (-dodge, sensitivity, "sensitivity_analysis", "10 features (cycling_fraction excluded)"),
    ]:
        est = data.loc[components, "proportion"].to_numpy()
        ci_low = data.loc[components, "ci_low"].to_numpy()
        ci_high = data.loc[components, "ci_high"].to_numpy()
        color = COLORS[key]
        ax.errorbar(
            est,
            y + offset,
            xerr=[est - ci_low, ci_high - est],
            fmt="o",
            color=color,
            markersize=8,
            markeredgewidth=0,
            elinewidth=2.0,
            capsize=0,
            label=label,
            zorder=2,
        )

    # Connecting segment between the two conditions' point estimates for
    # each component, making the shift itself the visual focus.
    for i, comp in enumerate(components):
        ax.plot(
            [original.loc[comp, "proportion"], sensitivity.loc[comp, "proportion"]],
            [y[i] + dodge, y[i] - dodge],
            color=COLORS["reference_line"],
            linewidth=1.3,
            zorder=0,
        )

    ax.set_yticks(y)
    ax.set_yticklabels([component_labels[c] for c in components])
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(len(components) - 0.5, -0.5)
    ax.set_xlabel("Proportion of variance (95% bootstrap CI)")
    ax.tick_params(axis="both", which="major", pad=5)
    # Legend placed below the axis, not inline: in this taller, narrower
    # canvas the CI lines span most of the plot area at every category,
    # leaving no empty quadrant for an inline legend without covering data.
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=1,
        bbox_to_anchor=(0.5, -0.06),
        handlelength=1.6,
        handletextpad=0.6,
    )
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    fig.savefig(output_path.with_suffix(".png"), dpi=600)
    plt.close(fig)
    return {"output_path": str(output_path)}


def _forest_panel(ax, labels: list[str], est, ci_low, ci_high, colors, *, xlabel: str) -> None:
    """Shared horizontal point-range ("forest plot") convention, used by
    every coefficient-style panel in Figure 5 so effect estimates are
    drawn identically throughout the manuscript. The CI here is primary
    content (not secondary context as in Figure 1's calibration curves),
    so it stays fully opaque and is sized for legibility, not de-emphasised."""
    y = np.arange(len(labels))
    est = np.asarray(est)
    ax.errorbar(
        est,
        y,
        xerr=[est - np.asarray(ci_low), np.asarray(ci_high) - est],
        fmt="none",
        ecolor=colors,
        elinewidth=2.1,
        capsize=0,
        zorder=1,
    )
    ax.scatter(est, y, c=colors, s=70, zorder=2, edgecolors="none")
    ax.axvline(0, color=COLORS["reference_line"], linewidth=1.1, linestyle="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="both", which="major", pad=5)


def build_figure_barrier_topology_with_ablation(project_root: Path) -> dict:
    model_path = project_root / "data" / "derived" / "barrier_topology_model_results.parquet"
    ablation_path = project_root / "data" / "derived" / "barrier_covariate_ablation.parquet"
    literature_path = project_root / "data" / "derived" / "literature_benchmark_results.parquet"
    output_path = (
        project_root / "reports" / "manuscript_figures" / "barrier_topology_with_ablation.pdf"
    )

    for p in (model_path, ablation_path, literature_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    model = pd.read_parquet(model_path).set_index("covariate")
    ablation = pd.read_parquet(ablation_path).set_index("step")
    literature = pd.read_parquet(literature_path).set_index("covariate")

    apply_publication_style()
    # 2x2 grid, not a single row of 3 (too wide once shrunk to fit a
    # word-processor page) or a column of 3 (too tall/long once
    # rendered at a legible font size). A: top-left, B: top-right,
    # C: bottom-left, with C's own legend moved into the otherwise
    # empty bottom-right cell. Top row is taller than the bottom row
    # since B has 5 categories vs. 2 each for A and C.
    # constrained_layout, not tight_layout(): tight_layout() ignores an
    # explicit GridSpec's own spacing and can silently re-compress rows
    # into each other; constrained_layout is GridSpec-aware and resolves
    # overlaps itself.
    fig = plt.figure(figsize=(13.4, 11.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], hspace=0.2, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_legend = fig.add_subplot(gs[1, 1])
    ax_legend.axis("off")

    # A: the two barrier covariates in the full, adjusted model.
    ax = ax_a
    covs = ["fibroblast_barrier_fraction", "suppressive_myeloid_barrier_fraction"]
    labels = ["Fibroblast", "Suppressive\nmyeloid"]
    colors = [COLORS["barrier_fibroblast"], COLORS["barrier_suppressive_myeloid"]]
    _forest_panel(
        ax,
        labels,
        model.loc[covs, "estimate"],
        model.loc[covs, "ci_low"],
        model.loc[covs, "ci_high"],
        colors,
        xlabel="Fixed-effect estimate\n(95% bootstrap CI)",
    )
    panel_title(ax, "Adjusted for state + niche")
    panel_label(ax, "A")

    # B: covariate-ablation stress test on the suppressive-myeloid effect --
    # the three block-level conditions Results reports, the single largest
    # individual covariate, and the full joint-adjustment model.
    ax = ax_b
    steps = [
        "barrier_only",
        "state_block_only",
        "niche_block_only",
        "+ niche_archetype_4_fraction",
        "full_state_and_niche",
    ]
    labels = [
        "Barrier only",
        "+ state composition",
        "+ niche composition",
        "+ niche archetype 4\n(largest single covariate)",
        "Full model\n(state + niche)",
    ]
    colors = [COLORS["not_significant"]] * 4 + [COLORS["barrier_suppressive_myeloid"]]
    step_data = ablation.loc[steps]
    _forest_panel(
        ax,
        labels,
        step_data["estimate"],
        step_data["ci_low"],
        step_data["ci_high"],
        colors,
        xlabel="Suppressive-myeloid estimate\n(95% CI)",
    )
    panel_title(ax, "Covariate-ablation stress test")
    panel_label(ax, "B")

    # C: this project's raw (unadjusted) correlations against the one
    # available published cross-cancer comparator. Legend moved into
    # the blank bottom-right grid cell rather than squeezed below C.
    ax = ax_c
    covs = ["suppressive_myeloid_barrier_fraction", "fibroblast_barrier_fraction"]
    labels = ["Suppressive\nmyeloid", "Fibroblast"]
    y = np.arange(len(covs))
    this_r = literature.loc[covs, "project_raw_r"].to_numpy()
    ax.scatter(
        this_r,
        y,
        c=COLORS["this_project"],
        s=70,
        zorder=2,
        edgecolors="none",
        label="This study (HNSCC, raw r)",
    )
    published_r = literature.loc["suppressive_myeloid_barrier_fraction", "published_r"]
    ax.scatter(
        [published_r],
        [0],
        c=COLORS["published_comparator"],
        s=70,
        marker="D",
        zorder=2,
        edgecolors="none",
        label="Grout et al. 2022 (lung, raw r)",
    )
    ax.axvline(0, color=COLORS["reference_line"], linewidth=1.1, linestyle="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_ylim(len(covs) - 0.5, -0.5)
    ax.set_xlabel("Raw (unadjusted) correlation")
    ax.tick_params(axis="both", which="major", pad=5)
    # panel_title()'s default 10pt pad sits too close to the axes once
    # the axes is this short (row 2 of the grid): the panel letter uses
    # axes-fraction positioning (a fixed 5% above the axes top), which
    # for a short axes is fewer absolute points than the title needs,
    # so the two collide. Setting the title directly with a larger
    # points-based pad gives the letter (still at its normal axes-
    # fraction offset) a clear gap to sit in below the title.
    ax.set_title(
        "Raw correlation vs. published estimate",
        fontsize=FS_PANEL_TITLE,
        fontweight="regular",
        pad=32,
        loc="center",
    )
    panel_label(ax, "C")
    handles, labels_ = ax.get_legend_handles_labels()
    ax_legend.legend(
        handles, labels_, loc="center", handlelength=1.6, fontsize=FS_LEGEND, frameon=False
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    fig.savefig(output_path.with_suffix(".png"), dpi=600)
    plt.close(fig)
    return {"output_path": str(output_path)}


def _lollipop_panel(ax, labels, values, *, threshold: float, xlabel: str) -> None:
    """Shared horizontal lollipop convention for ranked q-value panels --
    less ink than a filled bar for the same information, and visually
    consistent with the point-range idiom used throughout Figures 4-5."""
    y = np.arange(len(labels))
    ax.hlines(y, 0, values, color=COLORS["not_significant"], linewidth=1.6, zorder=1)
    ax.scatter(values, y, color=COLORS["not_significant"], s=42, zorder=2, edgecolors="none")
    ax.axvline(threshold, color=COLORS["reference_line"], linewidth=1.1, linestyle="--", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="both", which="major", pad=5)


def build_figure_hpv_consolidated(project_root: Path) -> dict:
    composition_path = (
        project_root / "data" / "derived" / "hpv_composition_comparison_results.parquet"
    )
    structure_path = project_root / "data" / "derived" / "hpv_structure_comparison_results.parquet"
    hpv_metadata_path = project_root / "metadata" / "hpv_status_validated.tsv"
    output_path = (
        project_root
        / "reports"
        / "manuscript_figures"
        / "hpv_composition_structure_and_discordance_qc.pdf"
    )

    for p in (composition_path, structure_path, hpv_metadata_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    composition = pd.read_parquet(composition_path).sort_values("pvalue_bh")
    structure = pd.read_parquet(structure_path).sort_values("pvalue_bh")
    hpv_status = pd.read_csv(hpv_metadata_path, sep="\t")

    apply_publication_style()
    # B (the longest ranked list, 13 categories) occupies its own full-
    # height column on the right, kept as a single unbroken list rather
    # than split into two half-columns -- a split list reads as two
    # panels sharing one letter, which is more confusing than a bit of
    # extra height. A and C (lighter content: 12 lineages, 6 patients)
    # stack top-to-bottom in the left column instead, at a height split
    # (12:6 rows) proportional to how much each actually needs.
    # constrained_layout, not tight_layout(): tight_layout() ignores an
    # explicit GridSpec's own spacing and can silently re-compress rows
    # into each other; constrained_layout is GridSpec-aware and resolves
    # overlaps itself.
    fig = plt.figure(figsize=(14.0, 13.0), constrained_layout=True)
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.75, 1], width_ratios=[1, 1.15], hspace=0.3, wspace=0.35
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_b = fig.add_subplot(gs[:, 1])
    threshold = -np.log10(0.05)

    lineage_display = {
        "Mast_cell": "Mast cell",
        "NK_cell": "NK cell",
        "Dendritic_cell": "Dendritic cell",
        "Perivascular_SmoothMuscle": "Perivascular / smooth muscle",
        "Epithelial_Tumour": "Epithelial (tumour)",
        "B_cell": "B cell",
        "Plasma_cell": "Plasma cell",
        "T_cell": "T cell",
    }
    ax = ax_a
    comp_labels = composition["lineage"].replace(lineage_display)
    _lollipop_panel(
        ax,
        comp_labels,
        -np.log10(composition["pvalue_bh"]),
        threshold=threshold,
        xlabel="-log10(BH q-value)",
    )
    # Direct ax.set_title() with an enlarged pad, not panel_title()'s
    # default pad=10: at this panel's row height, panel_label()'s
    # axes-fraction y-offset for the bold letter translates to fewer
    # absolute points than the default title padding needs, so the two
    # collide unless the title is pushed further up explicitly.
    ax.set_title(
        f"Cellular composition (n = {len(composition)})",
        fontsize=FS_PANEL_TITLE,
        fontweight="regular",
        pad=30,
        loc="center",
    )
    # Extra headroom above the topmost row ("Erythroid"): _lollipop_panel's
    # default half-row top pad isn't enough vertical room for the bold
    # panel letter above it, so the letter's own glyph was touching the
    # first row's label.
    ax.set_ylim(top=-1.2)
    # Default panel-letter x offset (-0.16) sits too close to this
    # narrower left-column panel's own y-tick category labels (e.g.
    # "Erythroid"); shifted further left to clear them.
    panel_label(ax, "A", x=-0.34)

    # Wrapped onto two lines: these category names are long enough that
    # a one-line label would force the y-axis tick-label margin to
    # consume most of the panel width, in turn crowding the panel
    # letter and title into the plot area itself.
    struct_labels = np.where(
        structure["category"] == structure["metric"],
        "Clone ecological-structure score",
        structure["category"] + "\n(" + structure["metric"].str.replace("_", " ") + ")",
    )
    ax = ax_b
    _lollipop_panel(
        ax,
        struct_labels,
        -np.log10(structure["pvalue_bh"]),
        threshold=threshold,
        xlabel="-log10(BH q-value)",
    )
    # Larger pad than A/C: panel_label()'s letter sits at a fixed
    # AXES-FRACTION offset above the axes top, which for this
    # full-height spanning column translates to more absolute points
    # than A/C's shorter axes -- enough to land inside the title's own
    # text height at the pad A/C use. A bigger pad keeps the letter
    # below the title's bottom edge, in the gap, as intended.
    ax.set_title(
        f"Ecosystem / clone structure (n = {len(structure)})",
        fontsize=FS_PANEL_TITLE,
        fontweight="regular",
        pad=50,
        loc="center",
    )
    panel_label(ax, "B")

    # C: per-patient clinical p16 vs. Xenium HPV16 E6/E7 probe signal --
    # an actual per-patient plot of the same comparison the previous
    # design summarised as bullet-point text.
    ax = ax_c
    tested = hpv_status.dropna(subset=["hpv_e6_e7_probe_positive_fraction"]).copy()
    tested = tested.sort_values("hpv_e6_e7_probe_positive_fraction", ascending=False)
    status_colors = {
        "confirmed_positive": COLORS["hpv_confirmed_positive"],
        "discordant_clinical_positive_probe_negative": COLORS["hpv_discordant"],
        "probe_positive_clinically_untested": COLORS["hpv_probe_positive_untested"],
    }
    status_display = {
        "confirmed_positive": "Confirmed positive",
        "discordant_clinical_positive_probe_negative": "Clinically positive,\nprobe-negative",
        "probe_positive_clinically_untested": "Probe-positive,\nclinically untested",
    }
    # Log-scale: a dot plot, not a lollipop -- "distance from zero" is not
    # a meaningful stem length once the x-axis is logarithmic.
    y = np.arange(len(tested))
    colors = tested["validated_hpv_status"].map(status_colors)
    ax.scatter(
        tested["hpv_e6_e7_probe_positive_fraction"], y, c=colors, s=64, zorder=2, edgecolors="none"
    )
    ax.set_yticks(y)
    ax.set_yticklabels(tested["patient_id"])
    ax.invert_yaxis()
    ax.set_ylim(len(tested) - 0.5, -0.5)
    # Extra headroom above the topmost row ("P09"): the default half-row
    # top pad isn't enough vertical room for the bold panel letter above
    # it, so the letter's own glyph was touching the first row's label.
    ax.set_ylim(top=-1.0)
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 2)
    ax.set_xlabel("HPV16 E6/E7 probe-\npositive fraction of cells")
    ax.tick_params(axis="both", which="major", pad=5)
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=c,
            markersize=11,
            label=status_display[k],
        )
        for k, c in status_colors.items()
    ]
    # Legend placed below the axis, not inside it: C is narrower and
    # shorter now that it's stacked under A rather than sharing a wide
    # top row, so no inside corner stays clear of the data across
    # its full x-range (probe signal and clinical status track
    # together here) -- both "lower right" and "upper left" ended up
    # sitting on top of points at this panel's new proportions.
    # bbox_to_anchor y pushed well past the two-line x-axis label
    # (-0.32 wasn't enough clearance and collided with it directly).
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.62),
        fontsize=FS_LEGEND - 1,
        handletextpad=0.5,
        framealpha=0.9,
        frameon=False,
    )
    ax.set_title(
        "Clinical p16 vs. probe signal",
        fontsize=FS_PANEL_TITLE,
        fontweight="regular",
        pad=30,
        loc="center",
    )
    # Default panel-letter x offset (-0.16) sits too close to this
    # narrower left-column panel's own y-tick patient labels (e.g.
    # "P09"); shifted further left to clear them.
    panel_label(ax, "C", x=-0.34)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    fig.savefig(output_path.with_suffix(".png"), dpi=600)
    plt.close(fig)
    return {"output_path": str(output_path)}
