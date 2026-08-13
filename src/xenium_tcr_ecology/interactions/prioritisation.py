"""Ranks `14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`'s candidate spatial interactions for
follow-up testing (`14_spatial_interactions_and_barriers/05_prioritise_testable_interactions.py`).

Exploratory synthesis, not a new statistical test -- combines results
already computed by earlier Spatial Interactions and Barriers milestones: spatial effect
size and cross-section/cross-patient consistency (`14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`'s
`spatial_interaction_scores.parquet`), panel completeness (Phase
14.01's `panel_supported_interactions.tsv`), and a qualitative
external-support tier (this module). This ranks candidate spatial
associations for prioritised follow-up -- it does not claim or imply
which, if any, causal pathway underlies an association; no claim of
this kind is made anywhere in this module's outputs.

**"External support" tier:** not a systematic literature database
query -- a qualitative tier grounded in already-established, well-known
immuno-oncology knowledge (the same knowledge already used to justify
each pair's inclusion in Phase 14.01's hand-curated list), stated
explicitly as qualitative rather than presented as a precise, sourced
literature count. `extensive`: axes with approved clinical therapeutics
and extensive solid-tumour spatial profiling (PD-1/PD-L1, CTLA-4,
foundational CD28 costimulation). `well_established`: mature
immuno-oncology axes with a shorter or more recent validation history
(TIGIT, LAG-3 -- relatlimab approved 2022; CCR7/CCL19 homing biology,
well-established but here applied outside its classic lymph-node
context). `moderate`: known but less extensively spatially profiled in
solid tumours (CD27/CD70).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

SIGNIFICANCE_THRESHOLD = 0.05

EXTERNAL_SUPPORT_TIERS: dict[str, str] = {
    "PD1_PDL1": "extensive",
    "CTLA4_CD80": "extensive",
    "CTLA4_CD86": "extensive",
    "CD28_CD80": "extensive",
    "CD28_CD86": "extensive",
    "TIGIT_PVR": "well_established",
    "LAG3_MHCII": "well_established",
    "CCR7_CCL19": "well_established",
    "CD27_CD70": "moderate",
}
EXTERNAL_SUPPORT_WEIGHTS: dict[str, float] = {
    "extensive": 1.0,
    "well_established": 0.75,
    "moderate": 0.5,
    "limited": 0.25,
}
COMPOSITE_CRITERIA = [
    "effect_size",
    "spatial_specificity",
    "cross_patient_consistency",
    "panel_completeness",
    "external_support_weight",
]


def aggregate_spatial_scores(
    scores: pd.DataFrame, section_to_patient: dict[str, str]
) -> pd.DataFrame:
    """Pure, testable: per-(sender_receiver_pair_id, lr_pair_id)
    aggregation of `14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`'s per-section results -- spatial
    specificity (fraction of tested sections significant),
    cross-patient consistency (fraction of distinct patients with at
    least one significant section, deduplicating replicate sections from
    the same patient so a 2-section patient cannot count twice), and an
    effect size (median log2 observed/null ratio across sections with
    both a positive observed score and a positive null mean -- rows
    where either is exactly zero are excluded from this ratio, not
    divided through with an arbitrary pseudocount)."""
    working = scores.copy()
    working["patient_id"] = working["section_id"].map(section_to_patient)
    working["is_significant"] = working["pvalue"] < SIGNIFICANCE_THRESHOLD

    rows = []
    for (sr_id, lr_id), group in working.groupby(["sender_receiver_pair_id", "lr_pair_id"]):
        n_sections_tested = len(group)
        n_sections_significant = int(group["is_significant"].sum())
        n_patients_tested = group["patient_id"].nunique()
        n_patients_with_significant_section = group.loc[
            group["is_significant"], "patient_id"
        ].nunique()

        valid_ratio_rows = group[(group["observed_mean_score"] > 0) & (group["null_mean"] > 0)]
        if len(valid_ratio_rows) > 0:
            log2_ratios = np.log2(
                valid_ratio_rows["observed_mean_score"] / valid_ratio_rows["null_mean"]
            )
            effect_size = float(log2_ratios.median())
        else:
            effect_size = float("nan")

        rows.append(
            {
                "sender_receiver_pair_id": sr_id,
                "lr_pair_id": lr_id,
                "n_sections_tested": n_sections_tested,
                "n_sections_significant": n_sections_significant,
                "spatial_specificity": n_sections_significant / n_sections_tested,
                "n_patients_tested": n_patients_tested,
                "n_patients_with_significant_section": n_patients_with_significant_section,
                "cross_patient_consistency": n_patients_with_significant_section
                / n_patients_tested,
                "effect_size": effect_size,
                "n_sections_effect_size": len(valid_ratio_rows),
            }
        )
    return pd.DataFrame(rows)


def compute_composite_priority(table: pd.DataFrame, criteria: list[str]) -> pd.DataFrame:
    """Pure, testable: transparent equal-weight rank-average composite
    across `criteria` (each ranked so higher values rank first), not
    presented as an optimised or validated weighting scheme -- a simple
    aggregation for prioritisation, not a new inferential statistic."""
    result = table.copy()
    percentile_cols = []
    for criterion in criteria:
        col = f"_{criterion}_percentile"
        result[col] = result[criterion].rank(pct=True, na_option="bottom")
        percentile_cols.append(col)
    result["composite_score"] = result[percentile_cols].mean(axis=1)
    result = result.drop(columns=percentile_cols)
    result = result.sort_values("composite_score", ascending=False).reset_index(drop=True)
    result["priority_rank"] = np.arange(1, len(result) + 1)
    return result


def build_interaction_priority_table(project_root: Path) -> dict:
    scores_path = project_root / "data" / "derived" / "spatial_interaction_scores.parquet"
    lr_pairs_path = project_root / "references" / "panel_supported_interactions.tsv"
    sample_manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    output_path = project_root / "results" / "interaction_priority_table.tsv"

    for path, phase in [
        (
            scores_path,
            "`14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`",
        ),
        (
            lr_pairs_path,
            "`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`",
        ),
        (sample_manifest_path, None),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found." + (f" Run {phase} first." if phase else ""))

    scores = pd.read_parquet(scores_path)
    lr_pairs = pd.read_csv(lr_pairs_path, sep="\t")
    sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
    section_to_patient = dict(zip(sample_manifest["section_id"], sample_manifest["patient_id"]))

    aggregated = aggregate_spatial_scores(scores, section_to_patient)
    aggregated = aggregated.merge(
        lr_pairs[["pair_id", "ligand", "receptor", "pair_complete"]],
        left_on="lr_pair_id",
        right_on="pair_id",
        how="left",
    )
    aggregated["panel_completeness"] = aggregated["pair_complete"].astype(float)

    missing_tier = sorted(set(aggregated["lr_pair_id"]) - set(EXTERNAL_SUPPORT_TIERS))
    if missing_tier:
        raise PipelineError(
            f"No external-support tier assigned for LR pair(s) {missing_tier} in EXTERNAL_SUPPORT_TIERS."
        )
    aggregated["external_support_tier"] = aggregated["lr_pair_id"].map(EXTERNAL_SUPPORT_TIERS)
    aggregated["external_support_weight"] = aggregated["external_support_tier"].map(
        EXTERNAL_SUPPORT_WEIGHTS
    )

    result = compute_composite_priority(aggregated, COMPOSITE_CRITERIA)
    result = result.drop(columns=["pair_id"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "n_candidate_interactions": len(result),
        "top_priority_sender_receiver_pair_id": str(result.iloc[0]["sender_receiver_pair_id"]),
        "top_priority_lr_pair_id": str(result.iloc[0]["lr_pair_id"]),
        "output_path": str(output_path),
    }
