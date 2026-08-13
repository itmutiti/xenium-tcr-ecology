"""Shared publication style for every main-text figure.

One typography, colour, and layout convention used across all six main
figures (`release/redesigned_main_figures.py`,
`tcr/vdj_ground_truth_validation.py`), so the manuscript reads as one
paper rather than six independently styled plots. Figure legends (the
interpretive text) live in the manuscript text, not inside the image:
figures here carry only a short, neutral panel descriptor -- never a
conclusion, decision, or verdict -- matching Genome Biology / Nature
Methods convention. This module only changes how already-computed
results are drawn; it does not read, compute, or alter any statistic.

The colour palette is Okabe & Ito (2008), the standard colourblind-safe
qualitative palette used throughout the Nature/Science-family journals.
Roles are assigned once, here, and reused by name everywhere so the same
concept (e.g. "spatial context", "suppressive-myeloid barrier") is
always the same colour across every figure it appears in.
"""

from __future__ import annotations

import matplotlib as mpl

# -- Okabe-Ito colourblind-safe palette (hex), named by role ---------------
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "grey": "#999999",
}

# Semantic colour roles, reused verbatim across every figure that touches
# these concepts, so a colour always means the same thing manuscript-wide.
COLORS = {
    # Three null-model constructions (Figure 1)
    "null_constrained_permutation": OKABE_ITO["blue"],
    "null_degree_preserving": OKABE_ITO["bluish_green"],
    "null_graph_preserving": OKABE_ITO["vermillion"],
    # Variance-partition components (Figure 4)
    "component_context": OKABE_ITO["bluish_green"],
    "component_patient": OKABE_ITO["blue"],
    "component_identity": OKABE_ITO["orange"],
    # Barrier topology (Figure 5)
    "barrier_suppressive_myeloid": OKABE_ITO["vermillion"],
    "barrier_fibroblast": OKABE_ITO["grey"],
    "this_project": OKABE_ITO["blue"],
    "published_comparator": OKABE_ITO["grey"],
    # TCR probe validation (Figure 2)
    "confirmed": OKABE_ITO["blue"],
    "not_confirmed": OKABE_ITO["grey"],
    # HPV status (Figure 6): per-patient clinical/molecular concordance
    "hpv_confirmed_positive": OKABE_ITO["blue"],
    "hpv_discordant": OKABE_ITO["vermillion"],
    "hpv_probe_positive_untested": OKABE_ITO["reddish_purple"],
    "hpv_test_value": OKABE_ITO["sky_blue"],
    # Generic
    "reference_line": "#B0B0B0",
    "not_significant": OKABE_ITO["grey"],
    "significant": OKABE_ITO["black"],
    # Primary analysis vs. a sensitivity/robustness variant thereof --
    # reused wherever a figure compares a prespecified analysis against a
    # named alternative (e.g. Figure 3's selected-K highlight, Figure 4's
    # feature-set sensitivity check).
    "primary_analysis": OKABE_ITO["blue"],
    "sensitivity_analysis": OKABE_ITO["vermillion"],
}

FONT_FAMILY = ["Arial", "Liberation Sans", "DejaVu Sans"]

# Typography scale (points). Revised 2026-07-19 (twice): first to a
# larger, more print-legible scale after reviewer feedback that the
# original scale (panel label 13 / axis label 9.5 / tick label 8.5 /
# title 9 / legend 8) read as undersized in print; then increased again
# after feedback that text was still too small once a wide multi-panel
# figure was pasted into a word processor and auto-scaled down to fit a
# page -- that scaling shrinks text in proportion to the figure's
# physical width, so legibility after a large scale-down depends on the
# font-to-canvas-width RATIO, not the absolute point size or export DPI.
# This second pass raises font size well ahead of canvas size to widen
# that ratio; canvas (`figsize`) values were deliberately left as-is.
FS_PANEL_LABEL = 26  # bold "A", "B", "C" panel letters
FS_AXIS_LABEL = 19
FS_TICK_LABEL = 17
FS_PANEL_TITLE = 19.5  # short, neutral, one line -- not a conclusion
FS_LEGEND = 17
FS_ANNOTATION = 16


def apply_publication_style() -> None:
    """Set matplotlib rcParams once, globally, before any figure is built.

    Idempotent -- safe to call at the top of every figure-building
    function without side effects across calls.
    """
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": FONT_FAMILY,
            "font.size": FS_TICK_LABEL,
            "axes.titlesize": FS_PANEL_TITLE,
            "axes.titleweight": "regular",
            "axes.labelsize": FS_AXIS_LABEL,
            "xtick.labelsize": FS_TICK_LABEL,
            "ytick.labelsize": FS_TICK_LABEL,
            "legend.fontsize": FS_LEGEND,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#4D4D4D",
            "axes.linewidth": 1.2,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.color": "#4D4D4D",
            "ytick.color": "#4D4D4D",
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "axes.grid": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,  # embed as editable, searchable text, not curves
            "figure.dpi": 300,
            "savefig.dpi": 600,  # print-reproduction minimum; PDF/SVG are vector and unaffected
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )


def panel_label(ax, letter: str, *, x: float = -0.16, y: float = 1.05) -> None:
    """Bold panel letter (A, B, C...) at a uniform position/size across
    every figure -- the only "label" a panel needs; the descriptive text
    belongs in the external figure legend, not repeated on the panel."""
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=FS_PANEL_LABEL,
        fontweight="bold",
        va="top",
        ha="left",
    )


def panel_title(ax, text: str) -> None:
    """Short, neutral, single-line panel descriptor (identifies the
    analysis; never states its result). Distinct from `panel_label`,
    which is the bold letter."""
    ax.set_title(text, fontsize=FS_PANEL_TITLE, fontweight="regular", pad=10, loc="center")


def significance_reference_line(ax, y: float = 0.0, **kwargs) -> None:
    """Thin, unobtrusive reference line (e.g. null-effect or alpha
    threshold) -- consistent style everywhere it is used."""
    style = {"color": COLORS["reference_line"], "linewidth": 0.9, "linestyle": "--", "zorder": 0}
    style.update(kwargs)
    ax.axhline(y, **style)
