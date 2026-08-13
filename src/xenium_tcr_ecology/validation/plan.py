"""Predeclares this project's external-validation claims, datasets and
success criteria before any External Validation and Generalisation dataset is acquired
or analysed (`16_external_validation_and_generalisation/00_define_validation_claims.py`).

**Sequencing requirement:** the registered `q1_framework_generalisation`
analysis (`governance/analysis_registry.tsv`) states its
exclusion/dataset criteria are "recorded in governance/validation_
plan.tsv (`16_external_validation_and_generalisation/00_define_validation_claims.py`) before this analysis runs" -- this module is
written and its output committed before `16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py` (dataset
acquisition) or `16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py` (the confirmatory test
itself) run, so the dataset-selection criteria below are predeclared,
not written retroactively once a specific dataset is already in hand.

**Claim set, each grounded in an already-established finding or
already-registered analysis:** `q1_framework_generalisation` reuses the
exact registered analysis_id and success criterion already defined in
`governance/analysis_registry.tsv` (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s calibration bounds, not a new threshold
chosen here). `cell_state_signature_generalisation` extends
`12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`-12.01's already-established methodology (module
coherence, state-proportion comparison against an independent
scRNA-seq reference) to whichever new reference(s) `16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`
acquires, using the same bar `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` already used (module
coherence significant in >=3/4 programmes) for direct comparability.
`ecosystem_signature_bulk_validation` and `source_paper_reproduction`
are new claims for External Validation and Generalisation specifically, each with an
explicitly cautious, qualitative-first success criterion (this
project's established discipline against overclaiming from
underpowered or indirect comparisons).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

VALIDATION_CLAIMS: list[dict] = [
    {
        "claim_id": "q1_framework_generalisation",
        "claim": "The calibrated spatial-ecology framework (null models, rarefied descriptors, graph pruning) produces consistent, correctly-calibrated results when applied end-to-end to an independent spatial dataset, not only this project's own GSE300147/McCord et al. cohort.",
        "validation_dataset": "An independent spatial transcriptomics dataset acquired in `16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py`, publicly available, with per-cell spatial coordinates and either clonal/lineage information or a proxy cell-type annotation sufficient to re-run the null-model calibration suite.",
        "validation_method": "Re-run `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s already-calibrated null-model calibration suite (constrained permutation, degree-preserving, graph-preserving) end-to-end on the new dataset (`16_external_validation_and_generalisation/05_validate_framework_on_independent_dataset.py`).",
        "success_criterion": "Empirical Type I error at nominal alpha=0.05 falls within the same bounds `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py` established on this project's synthetic ground truth -- reusing the exact registered analysis_id q1_framework_generalisation and its stated success criterion in governance/analysis_registry.tsv, not a new threshold invented here.",
        "phase_reference": "16.01, 16.05",
    },
    {
        "claim_id": "cell_state_signature_generalisation",
        "claim": "This project's T-cell state and checkpoint/chemokine/interferon/antigen-presentation programme gene signatures (`05_preprocessing_and_normalisation/03_calculate_program_scores.py`, Cell Type Annotation) correspond to coherent states in an independent full-transcriptome HNSCC reference, not only the targeted 623-gene Xenium panel.",
        "validation_dataset": "Public HNSCC scRNA-seq reference(s) acquired in `16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`, full-transcriptome, with T-cell/stromal annotations.",
        "validation_method": "Reuse `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`, `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`'s already-validated methodology (module coherence via mean pairwise correlation against a matched background gene pool, state-proportion comparison) on the new reference(s) (`16_external_validation_and_generalisation/03_validate_cell_state_signatures.R`).",
        "success_criterion": "Module coherence significant (p<0.05) in at least 3/4 of the tested programmes -- the same bar `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py` already used on GSE103322, for direct comparability across the two independent references.",
        "phase_reference": "16.02, 16.03",
    },
    {
        "claim_id": "ecosystem_signature_bulk_validation",
        "claim": "This project's Niche and Ecosystem Discovery/ecosystem definitions correspond to detectable signal in an independent bulk RNA-seq HNSCC cohort with clinical/immune annotation.",
        "validation_dataset": "A public bulk RNA-seq HNSCC cohort with clinical and/or immune-infiltration annotation (e.g. TCGA-HNSC or a comparable public cohort), acquired in `16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`.",
        "validation_method": "Project ecosystem-derived gene signatures onto the bulk cohort (e.g. ssGSEA-style scoring) and test correlation with clinical/immune annotation (`16_external_validation_and_generalisation/04_validate_ecosystem_signatures_in_bulk.py`).",
        "success_criterion": "Qualitative first: signature scores correlate in the biologically-expected direction with at least one independent immune-infiltration proxy already established as meaningful in the bulk-RNA-seq literature; any quantitative claim is reported with effect size and uncertainty, not treated as confirmatory (this outcome domain has no prespecified confirmatory status in governance/analysis_registry.tsv).",
        "phase_reference": "16.02, 16.04",
    },
    {
        "claim_id": "source_paper_reproduction",
        "claim": "Selected, previously-published McCord et al. findings are reproduced by this project's independent reanalysis pipeline.",
        "validation_dataset": "No new dataset -- this project's primary cohort, compared directly against the source paper's already-published statistics.",
        "validation_method": "Direct, claim-by-claim comparison against the source paper's already-published statistics (`16_external_validation_and_generalisation/06_compare_with_source_paper_results.py`).",
        "success_criterion": "Qualitative concordance in direction for each pre-selected claim (e.g. replicate reproducibility's 1-of-7-discordant-pairs finding, already independently reproduced in `04_quality_control/08_assess_replicate_concordance.R`); any quantitative discrepancy is reported with its magnitude.",
        "phase_reference": "16.06",
    },
]


def validate_claims_well_formed(claims: list[dict]) -> None:
    """Structural check: every claim has all required, non-empty
    fields, and every `claim_id` is unique."""
    required_fields = [
        "claim_id",
        "claim",
        "validation_dataset",
        "validation_method",
        "success_criterion",
        "phase_reference",
    ]
    seen_ids = set()
    for claim in claims:
        for field in required_fields:
            if not claim.get(field):
                raise PipelineError(
                    f"Claim '{claim.get('claim_id', '<unknown>')}' is missing a non-empty '{field}'."
                )
        if claim["claim_id"] in seen_ids:
            raise PipelineError(f"Duplicate claim_id '{claim['claim_id']}'.")
        seen_ids.add(claim["claim_id"])


def build_validation_plan(project_root: Path) -> dict:
    output_path = project_root / "governance" / "validation_plan.tsv"

    validate_claims_well_formed(VALIDATION_CLAIMS)

    result = pd.DataFrame(VALIDATION_CLAIMS)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "n_claims": len(result),
        "output_path": str(output_path),
    }
