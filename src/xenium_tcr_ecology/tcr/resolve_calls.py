"""Resolve multiclonal and ambiguous TCR calls (`08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`).

Classifies every T cell evaluated in `08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py` into exactly one of four
categories, combining evidence already built in this phase rather than
re-deriving it: `unassigned`, `singlet`, `probable_multiplet`,
`low_confidence`.

**Classification precedence, in order:**
1. `unassigned` -- no probe detected at all (`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s `any_detection
   == False`).
2. `probable_multiplet` -- ambiguous multi-probe detection
   (`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`'s `is_multi_probe_ambiguous_excluding_likely_pairs`, i.e.
   not explained by a normal single-clone TRA+TRB pair).
3. `low_confidence` -- the detected probe(s) include at least one with a
   high empirical false-positive-rate estimate (`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`,
   `empirical_fpr > LOW_CONFIDENCE_FPR_THRESHOLD`) -- a technically
   "clean" singlet call built on an unreliable probe should not be
   reported with the same confidence as one built on a well-validated
   probe.
4. `singlet` -- everything else: exactly one detected probe, or a
   TRA+TRB pair consistent with one normally-paired clone, using only
   probes with an acceptable empirical false-positive rate.

**Cross-check against an independent QC signal, not part of the
classification decision itself:** `probable_multiplet` cells are compared
against `04_quality_control/04_estimate_transcript_spillover.py`'s per-cell `spillover_risk_score` (cell-boundary-
polygon-based transcript-spillover risk, entirely independent of anything
computed in TCR Clonal Analysis) as an external validity check -- if
ambiguous multi-probe cells are enriched for high spillover risk relative
to singlets, that is independent corroboration that at least some
`probable_multiplet` calls reflect a genuine segmentation/spillover
artefact rather than two truly co-resident clones in one cell, without
claiming this classification scheme resolves that distinction (it does
not attempt to for individual cells).

`LOW_CONFIDENCE_FPR_THRESHOLD = 0.5`: checked against the
`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R` distribution before choosing -- 23/105 probes (21.9%) exceed
this, a meaningful minority, and the threshold itself has a direct
interpretation ("at least half of this probe's positive calls could
plausibly be background noise alone").
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

LOW_CONFIDENCE_FPR_THRESHOLD = 0.5


def classify_calls(
    calls: pd.DataFrame,
    probe_fpr: dict[str, float],
    fpr_threshold: float = LOW_CONFIDENCE_FPR_THRESHOLD,
) -> pd.Series:
    """Pure, testable classification. `calls` must have `any_detection`,
    `is_multi_probe_ambiguous_excluding_likely_pairs`, and
    `detected_probes` (semicolon-joined probe names) columns."""

    def _classify_row(row) -> str:
        if not row["any_detection"]:
            return "unassigned"
        if row["is_multi_probe_ambiguous_excluding_likely_pairs"]:
            return "probable_multiplet"
        probes = [p for p in row["detected_probes"].split(";") if p]
        max_fpr = max((probe_fpr.get(p, 0.0) for p in probes), default=0.0)
        if max_fpr > fpr_threshold:
            return "low_confidence"
        return "singlet"

    return calls.apply(_classify_row, axis=1)


def build_resolved_calls(project_root: Path) -> dict:
    calls_path = project_root / "data" / "derived" / "tcr_cell_calls.parquet"
    fpr_path = project_root / "data" / "derived" / "tcr_false_positive_estimates.parquet"
    spillover_path = project_root / "data" / "derived" / "spillover_risk.parquet"
    output_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"

    for p in (calls_path, fpr_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    calls = pd.read_parquet(calls_path)
    fpr = pd.read_parquet(fpr_path)
    probe_fpr = fpr.set_index("probe_name")["empirical_fpr"].to_dict()

    resolved = calls.copy()
    resolved["resolution"] = classify_calls(calls, probe_fpr)

    spillover_comparison = None
    if spillover_path.is_file():
        spillover = pd.read_parquet(spillover_path)
        joined = resolved.join(spillover[["spillover_risk_score"]], how="inner")
        spillover_comparison = (
            joined.groupby("resolution", observed=True)["spillover_risk_score"]
            .mean()
            .round(4)
            .to_dict()
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    resolved.to_parquet(output_path)

    return {
        "n_cells": len(resolved),
        "resolution_counts": resolved["resolution"].value_counts().to_dict(),
        "fpr_threshold": LOW_CONFIDENCE_FPR_THRESHOLD,
        "mean_spillover_risk_by_resolution": spillover_comparison,
        "output_path": str(output_path),
    }
