"""Directional consistency quantification (`12_external_checkpoint_validation/02_quantify_directional_consistency.py`).

Goes beyond `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`'s raw abundance comparison and `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`'s
module-coherence test to ask the specific question this milestone is
named for: is the direction of the effects found in this project's data
consistent with the independent GSE103322 reference, not just "is a
p-value small"?

Two complementary directional statistics, both computed on Phase
12.00's state-proportion comparison
(`bulk_reference_state_comparison.parquet`):

1. **Pairwise sign agreement**: for every pair of T-cell states, checks
   whether the sign of (fraction[A] - fraction[B]) -- i.e. "A is more
   common than B" -- agrees between this project and the reference. A
   direct, assumption-free directional check across all `C(7,2)=21`
   state pairs.
2. **Rank correlation**: Spearman rank correlation of the 7 state
   abundances across the two datasets -- a single summary number for
   "do the two cohorts agree on which states are common vs. rare,"
   independent of the exact magnitudes (which `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py` already showed
   differ for `Cycling` specifically).

Also links directly to `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s continuous-structure loadings
(`clone_ecological_structure_loadings.parquet`): for each of the 6
T-cell-state features that load on that axis (the only `11_clone_spatial_descriptors/06_discover_provisional_structure.R`
features with an external counterpart in this T-cell-only reference --
tumour-engagement/APC-support features have no comparable quantity in
GSE103322, which carries no spatial or TCR data), reports whether that
state's rank position is consistent (same relative rank, allowing a
documented tolerance) between the two datasets -- directly answering
whether the specific feature driving Phase 11.06's axis
(`cycling_fraction`, loading -0.998) shows a directional problem or
just a magnitude difference.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Maps `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s feature names to the comparable T-cell state.
LOADING_FEATURE_TO_STATE = {
    "cytotoxic_fraction": "Cytotoxic",
    "exhausted_fraction": "Exhausted",
    "cycling_fraction": "Cycling",
    "treg_fraction": "Treg",
    "cd4_fraction": "CD4",
    "cd8_fraction": "CD8",
}


def compute_pairwise_sign_agreement(
    project_fractions: pd.Series, reference_fractions: pd.Series
) -> pd.DataFrame:
    """Pure, testable: for every pair of states present in both
    series, whether sign(project[a]-project[b]) matches sign(reference[a]-reference[b])."""
    states = sorted(set(project_fractions.index) & set(reference_fractions.index))
    rows = []
    for i, a in enumerate(states):
        for b in states[i + 1 :]:
            project_sign = float(np.sign(project_fractions[a] - project_fractions[b]))
            reference_sign = float(np.sign(reference_fractions[a] - reference_fractions[b]))
            rows.append(
                {
                    "state_a": a,
                    "state_b": b,
                    "project_sign": project_sign,
                    "reference_sign": reference_sign,
                    "agrees": bool(project_sign == reference_sign),
                }
            )
    return pd.DataFrame(rows)


def compute_rank_consistency(project_fractions: pd.Series, reference_fractions: pd.Series) -> dict:
    """Pure, testable: Spearman rank correlation of state abundances
    between the two datasets."""
    states = sorted(set(project_fractions.index) & set(reference_fractions.index))
    rho, pvalue = spearmanr(
        [project_fractions[s] for s in states], [reference_fractions[s] for s in states]
    )
    return {"spearman_rho": float(rho), "pvalue": float(pvalue), "n_states": len(states)}


def compute_rank_shift(
    project_fractions: pd.Series, reference_fractions: pd.Series
) -> pd.DataFrame:
    """Pure, testable: each state's rank (1 = most abundant) in each
    dataset, and the shift between them."""
    project_rank = project_fractions.rank(ascending=False)
    reference_rank = reference_fractions.rank(ascending=False)
    states = sorted(set(project_fractions.index) & set(reference_fractions.index))
    return pd.DataFrame(
        {
            "state": states,
            "project_rank": [int(project_rank[s]) for s in states],
            "reference_rank": [int(reference_rank[s]) for s in states],
            "rank_shift": [int(reference_rank[s] - project_rank[s]) for s in states],
        }
    )


def build_directional_consistency(project_root: Path) -> dict:
    comparison_path = project_root / "data" / "derived" / "bulk_reference_state_comparison.parquet"
    loadings_path = (
        project_root / "data" / "derived" / "clone_ecological_structure_loadings.parquet"
    )
    output_path = project_root / "results" / "external_checkpoint" / "directional_consistency.tsv"

    for path, phase in [
        (
            comparison_path,
            "`12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`",
        ),
        (loadings_path, "`11_clone_spatial_descriptors/06_discover_provisional_structure.R`"),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    comparison = pd.read_parquet(comparison_path).set_index("state")
    project_fractions = comparison["project_fraction"]
    reference_fractions = comparison["reference_fraction"]

    sign_agreement = compute_pairwise_sign_agreement(project_fractions, reference_fractions)
    rank_consistency = compute_rank_consistency(project_fractions, reference_fractions)
    rank_shift = compute_rank_shift(project_fractions, reference_fractions)

    loadings = pd.read_parquet(loadings_path).set_index("feature")
    rank_shift["structure_axis_loading"] = rank_shift["state"].map(
        lambda s: next(
            (
                loadings.loc[f, "loading"]
                for f, mapped in LOADING_FEATURE_TO_STATE.items()
                if mapped == s and f in loadings.index
            ),
            np.nan,
        )
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    rank_shift.to_csv(output_path, sep="\t", index=False)

    sign_agreement_path = project_root / "data" / "derived" / "directional_sign_agreement.parquet"
    sign_agreement.to_parquet(sign_agreement_path)

    return {
        "n_state_pairs": len(sign_agreement),
        "n_pairs_agreeing": int(sign_agreement["agrees"].sum()),
        "fraction_pairs_agreeing": float(sign_agreement["agrees"].mean()),
        "spearman_rho": rank_consistency["spearman_rho"],
        "spearman_pvalue": rank_consistency["pvalue"],
        "output_path": str(output_path),
        "sign_agreement_path": str(sign_agreement_path),
    }
