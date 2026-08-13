"""Freeze-or-revise decision (`12_external_checkpoint_validation/03_decide_freeze_or_revise.py`).

Applies a decision rule to `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`-12.02's external
sanity-check outputs. This rule was finalised in the same working
session as the checks it gates, rather than declared strictly ahead of
seeing their results; see `docs/analysis_amendments.md` (2026-07-11
entry) for that note. Every threshold used is a standard, generic,
field-convention value (`p < 0.05`; a positive sign for rank
correlation), not fitted to reproduce a specific outcome.

**Rule, two necessary conditions, both must hold for FREEZE:**
1. **Module transfer (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`):** all tested transcriptional
   programmes must clear their empirical null in both datasets
   (`transfers == True` for every row) -- if the underlying gene
   signature itself does not reproduce externally at all, no amount of
   directional nuance justifies freezing.
2. **Directional sign (`12_external_checkpoint_validation/02_quantify_directional_consistency.py`):** the overall Spearman rank
   correlation between project and reference state abundances must be
   positive -- the two cohorts must not be pointing in opposite
   directions overall. A `p < 0.05` bar is not required here: with only
   n=7 states, achieving conventional significance is statistically
   difficult regardless of true effect size, so requiring it would set an
   unreasonably strict bar unrelated to the question being asked (is the
   direction right at all, not "is it significant at this small n").

If both hold, the outcome is FREEZE, with an explicit caveat recorded
(not silently dropped) for any individual state showing a large rank
shift (`|rank_shift| >= REVISE_FLAG_RANK_SHIFT`) --
`12_external_checkpoint_validation/02_quantify_directional_consistency.py` already found `Cycling` and `Ambiguous` meet this, and that
caveat is carried forward into the freeze record rather than allowed to
vanish once the top-line decision is made.

REVISE (v2, capped at one revision) applies only if either necessary
condition fails.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

REVISE_FLAG_RANK_SHIFT = 3


def decide_freeze_or_revise(all_programs_transfer: bool, overall_spearman_rho: float) -> str:
    """Pure, testable: the two-condition decision rule."""
    if all_programs_transfer and overall_spearman_rho > 0:
        return "freeze"
    return "revise"


def identify_flagged_states(
    rank_shift_df: pd.DataFrame, threshold: int = REVISE_FLAG_RANK_SHIFT
) -> list[str]:
    """Pure, testable: states whose rank shift exceeds `threshold` in
    absolute value -- carried forward as an explicit caveat, not
    silently dropped once the top-line decision is made."""
    flagged = rank_shift_df[rank_shift_df["rank_shift"].abs() >= threshold]
    return sorted(flagged["state"].tolist())


def build_freeze_decision(project_root: Path) -> dict:
    program_transfer_path = project_root / "data" / "derived" / "program_transfer_results.parquet"
    directional_consistency_path = (
        project_root / "results" / "external_checkpoint" / "directional_consistency.tsv"
    )
    taxonomy_log_path = project_root / "governance" / "taxonomy_version_log.tsv"
    output_path = project_root / "governance" / "freeze_decision.tsv"

    for path, phase in [
        (
            program_transfer_path,
            "`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`",
        ),
        (
            directional_consistency_path,
            "`12_external_checkpoint_validation/02_quantify_directional_consistency.py`",
        ),
        (
            taxonomy_log_path,
            "`11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    program_transfer = pd.read_parquet(program_transfer_path)
    all_programs_transfer = bool(program_transfer["transfers"].all())

    rank_shift = pd.read_csv(directional_consistency_path, sep="\t")
    from scipy.stats import spearmanr

    rho, pvalue = spearmanr(rank_shift["project_rank"], rank_shift["reference_rank"])
    # project_rank/reference_rank are ranks (1=most abundant); correlating
    # ranks directly gives the same Spearman rho as correlating the
    # underlying fractions (`12_external_checkpoint_validation/02_quantify_directional_consistency.py`'s `compute_rank_consistency`),
    # reused here without re-reading the raw fraction data a second time.
    overall_spearman_rho = float(rho)

    decision = decide_freeze_or_revise(all_programs_transfer, overall_spearman_rho)
    flagged_states = identify_flagged_states(rank_shift)

    taxonomy_log = pd.read_csv(taxonomy_log_path, sep="\t")
    current_version = taxonomy_log.iloc[-1]["taxonomy_version"]
    next_version = current_version if decision == "freeze" else "v2_revised"

    entry = pd.DataFrame(
        [
            {
                "taxonomy_version_evaluated": current_version,
                "decision": decision,
                "resulting_version": next_version,
                "all_programs_transfer": all_programs_transfer,
                "overall_spearman_rho": overall_spearman_rho,
                "overall_spearman_pvalue": float(pvalue),
                "flagged_states": ";".join(flagged_states),
                "rule": "freeze iff (all `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` programs transfer) AND (`12_external_checkpoint_validation/02_quantify_directional_consistency.py` overall Spearman rho > 0); revise otherwise",
            }
        ]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    entry.to_csv(output_path, sep="\t", index=False)

    return {
        "decision": decision,
        "resulting_version": next_version,
        "all_programs_transfer": all_programs_transfer,
        "overall_spearman_rho": overall_spearman_rho,
        "flagged_states": flagged_states,
        "output_path": str(output_path),
    }
