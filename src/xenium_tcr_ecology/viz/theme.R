# Shared publication theme for every R-generated figure (main and
# supplementary), mirroring `xenium_tcr_ecology/viz/style.py`'s
# typography, colour-role, and layout conventions so R- and
# Python-generated figures read as one paper. Sourced, not a package --
# matches this codebase's existing convention of self-contained scripts
# (no R package/NAMESPACE exists here); every script that plots sources
# this one file by relative path from the project root.
#
# Only changes how already-computed results are drawn. Does not read,
# compute, or alter any statistic.

suppressPackageStartupMessages(library(ggplot2))

# Okabe & Ito (2008) colourblind-safe qualitative palette -- identical
# hex values to OKABE_ITO in style.py.
OKABE_ITO <- list(
  black = "#000000", orange = "#E69F00", sky_blue = "#56B4E9",
  bluish_green = "#009E73", yellow = "#F0E442", blue = "#0072B2",
  vermillion = "#D55E00", reddish_purple = "#CC79A7", grey = "#999999"
)

# Semantic colour roles reused across figures -- kept in sync with
# COLORS in style.py wherever the same concept appears in both R- and
# Python-generated figures.
PUB_COLORS <- list(
  primary_analysis = OKABE_ITO$blue,
  sensitivity_analysis = OKABE_ITO$vermillion,
  reference_line = "#B0B0B0",
  not_significant = OKABE_ITO$grey,
  significant = OKABE_ITO$black,
  confirmed = OKABE_ITO$blue,
  not_confirmed = OKABE_ITO$grey
)

# base_size raised twice (2026-07-19): 9.5->14 to match style.py's
# first FS_* increase, then 14->19 after feedback that text was still
# too small once a wide multi-panel figure was pasted into a word
# processor and auto-scaled to fit a page -- that shrinks text in
# proportion to the figure's physical width, so legibility depends on
# the font-to-canvas-width ratio, not absolute point size. Canvas
# (`open_publication_pdf` width/height) values were deliberately left
# as-is; only this ratio changed. All theme_minimal()-derived element
# sizes below cascade from base_size, so this one change rescales axis
# text, titles, strip text, and legend text together, proportionally.
theme_publication <- function(base_size = 19) {
  theme_minimal(base_family = "Liberation Sans", base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "grey30", linewidth = 0.6),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.6),
      axis.ticks.length = unit(4, "pt"),
      plot.title = element_text(size = base_size, face = "plain", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, face = "plain", colour = "grey30", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "plain"),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(12, 14, 8, 8)
    )
}

# Bold panel letter, positioned identically to panel_label() in
# style.py. Call after printing the plot to the current page/viewport,
# or pass letter= to annotate_panel_letter(plot, "A") to compose it into
# the ggplot object itself before printing.
annotate_panel_letter <- function(plot, letter) {
  plot + labs(tag = letter) +
    theme(
      plot.tag = element_text(size = 26, face = "bold", family = "Liberation Sans"),
      plot.tag.position = c(0.02, 0.98)
    )
}

# cairo_pdf(..., family="Liberation Sans") is required (not the base pdf()
# device) for the Liberation Sans family to actually render -- base pdf()
# only supports a handful of built-in PostScript font families.
open_publication_pdf <- function(path, width, height) {
  grDevices::cairo_pdf(path, width = width, height = height, family = "Liberation Sans", onefile = TRUE)
}

#' Lay out 1-6 ggplot objects on one page as lettered panels (A, B, C...),
#' in a single row or a wrapped grid, using base `grid` viewports (no
#' layout package -- patchwork/cowplot/gridExtra are not installed in
#' this environment). Replaces the previous per-script pattern of one
#' `print(p)` per page in a multi-page PDF: every multi-panel
#' supplementary figure is now one figure, not N stacked pages.
#' @param plots named or unnamed list of ggplot objects
#' @param ncol number of columns; defaults to one row (length(plots) columns)
compose_panels <- function(plots, ncol = NULL, widths = NULL, heights = NULL) {
  n <- length(plots)
  if (is.null(ncol)) ncol <- n
  nrow <- ceiling(n / ncol)
  if (is.null(widths)) widths <- rep(1, ncol)
  if (is.null(heights)) heights <- rep(1, nrow)

  grid::grid.newpage()
  layout <- grid::grid.layout(nrow = nrow, ncol = ncol, widths = unit(widths, "null"), heights = unit(heights, "null"))
  grid::pushViewport(grid::viewport(layout = layout))

  letters_seq <- LETTERS[seq_len(n)]
  for (i in seq_len(n)) {
    row <- ceiling(i / ncol)
    col <- ((i - 1) %% ncol) + 1
    grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
    print(plots[[i]], vp = grid::viewport())
    grid::grid.draw(grid::textGrob(
      letters_seq[i],
      x = unit(3, "pt"), y = unit(1, "npc") - unit(3, "pt"),
      just = c("left", "top"),
      gp = grid::gpar(fontface = "bold", fontsize = 26, fontfamily = "Liberation Sans")
    ))
    grid::popViewport()
  }
  grid::popViewport()
}
