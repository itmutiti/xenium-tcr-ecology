#!/usr/bin/env python3
"""
`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`

Tests whether each marker-gene programme (cytotoxicity, exhaustion,
proliferation, Treg) is a coherent co-expressed module in both
this project's T cells and the independent GSE103322 reference --
directly investigating `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`'s Cycling-state discrepancy
-- see src/xenium_tcr_ecology/external_checkpoint/program_transfer.py's module
docstring.

Primary output: reports/external_checkpoint/program_transfer.pdf
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
from xenium_tcr_ecology.external_checkpoint.program_transfer import build_program_transfer_test


def render_program_transfer_report(result: pd.DataFrame, output_path) -> None:
    from xenium_tcr_ecology.viz.style import COLORS, apply_publication_style, panel_title

    apply_publication_style()
    # Narrower, taller canvas: a word processor auto-scales a pasted
    # image to fit the page's text width (~6.5in), shrinking every
    # dimension including font size by that ratio. A narrower canvas
    # shrinks less, so the same in-image point size reads larger once
    # pasted.
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    y = np.arange(len(result))
    ax.hlines(
        y,
        result["project_observed_coherence"],
        result["reference_observed_coherence"],
        color=COLORS["reference_line"],
        linewidth=1.3,
        zorder=1,
    )
    ax.scatter(
        result["project_observed_coherence"],
        y,
        color=COLORS["primary_analysis"],
        s=60,
        zorder=2,
        edgecolors="none",
        label="This study (Xenium)",
    )
    ax.scatter(
        result["reference_observed_coherence"],
        y,
        color=COLORS["sensitivity_analysis"],
        s=60,
        zorder=2,
        marker="D",
        edgecolors="none",
        label="GSE103322 (external reference)",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(result["program"].str.capitalize())
    ax.axvline(0, color=COLORS["reference_line"], linewidth=1.1, linestyle="--", zorder=0)
    ax.set_xlabel("Mean pairwise Spearman correlation\namong programme marker genes")
    ax.tick_params(axis="both", which="major", pad=5)
    panel_title(ax, "Transcriptional programme coherence vs. external reference")
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
        script_name="01_test_transcriptional_program_transfer",
        project_root=project_root,
        phase="12_external_checkpoint_validation",
    )

    try:
        summary = build_program_transfer_test(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    result = pd.read_parquet(summary["output_path"])
    report_path = project_root / "reports" / "external_checkpoint" / "program_transfer.pdf"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    render_program_transfer_report(result, report_path)
    summary["report_path"] = str(report_path)

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_programs_tested']} programme(s) tested, {summary['n_programs_transferring']} transfer "
        f"(both p<0.05). Not transferring: {summary['programs_not_transferring']}. "
        f"Wrote {summary['output_path']}, {summary['report_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
