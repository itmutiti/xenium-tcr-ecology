"""Compares this project's completed findings against each of the 11
McCord et al. source-paper claims, and states which findings are new
(`16_external_validation_and_generalisation/06_compare_with_source_paper_results.py`) -- the `source_paper_reproduction`
claim (`16_external_validation_and_generalisation/00_define_validation_claims.py`).

**Direct claim-by-claim comparison:** all 11 rows below are addressed
explicitly, using already-completed results (not re-derived here) --
covering: `reproduced_and_extended` findings
(cohort-wide, null-model-calibrated, rarefaction-normalised quantitative
extensions of the source paper's single-patient/illustrative findings),
a `methodological_choice` (Niche and Ecosystem Discovery's
composition-vector consensus clustering vs. the source paper's
DBSCAN+Moran's I), one `diverges` finding under a literal same-metric
comparison (`04_quality_control/08_assess_replicate_concordance.R`'s
replicate concordance -- 0/7 pairs discordant under the paper's
CDR3-probe-Jaccard metric, vs. the paper's "1 of 7" finding), a substantive scope
divergence (`11_clone_spatial_descriptors/05_test_discrete_vs_
continuous_structure.R`'s prespecified continuous-structure finding vs.
Figure 2's 14 discrete transcriptional clusters), an
`addressed_limitation` (this project's HPV-Stratified Analysis/16 as a
direct response to the source paper's named n=10 caveat), and claims
`not_independently_tested` (Figures 1 and 3, which need flow-cytometry/
PBMC data this cohort's Xenium-only data does not have).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

VALID_STATUSES = {
    "reproduced",
    "reproduced_and_extended",
    "diverges",
    "methodological_choice",
    "addressed_limitation",
    "not_independently_tested",
}

SOURCE_PAPER_COMPARISON: list[dict] = [
    {
        "figure_or_table": "Figure 1",
        "comparison_status": "not_independently_tested",
        "this_project_reference": "n/a",
        "note": "Flow-cytometry-derived immune-composition groups and their p16/smoking/site clinical associations require flow cytometry data this cohort does not have.",
        "novel_contribution": "n/a -- explicitly out of scope.",
    },
    {
        "figure_or_table": "Figure 2",
        "comparison_status": "diverges",
        "this_project_reference": "`11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R` (scripts/11_clone_spatial_descriptors/05_test_discrete_vs_continuous_structure.R, prespecified)",
        "note": "The source paper finds 14 discrete transcriptional T-cell clusters (full-transcriptome scRNA-seq). This project's prespecified structure test (BIC + gap statistic, on spatial/clone-level descriptors from the targeted 623-gene panel) found continuous structure, not discrete clusters -- an unreconciled divergence.",
        "novel_contribution": "A quantitative, prespecified test of discrete-vs-continuous structure -- the source paper's 14-cluster claim was not independently statistically tested there, only presented as the output of standard scRNA clustering. A plausible reason for the divergence, not asserted as proven: a targeted 623-gene panel and spatial (not purely transcriptional) descriptors are a different measurement space from full-transcriptome scRNA clustering.",
    },
    {
        "figure_or_table": "Figure 3",
        "comparison_status": "not_independently_tested",
        "this_project_reference": "n/a",
        "note": "PBMC-vs-TIL clonal overlap (Morisita-Horn index) requires paired PBMC samples this cohort does not have. The bystander/viral-reactive-clone concept is acknowledged in this project's TCR Clonal Analysis clone ascertainment design, but not independently re-tested.",
        "novel_contribution": "n/a -- explicitly out of scope.",
    },
    {
        "figure_or_table": "Figure 4",
        "comparison_status": "reproduced_and_extended",
        "this_project_reference": "`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py` (tumour_engagement), 11.03 (apc_support/antigen_presentation), 14.03 (q3_barrier_topology_confirmatory, prespecified)",
        "note": "The source paper shows a qualitative, illustrative spatial exhaustion/naive gradient by proximity to tumour/myeloid/B-cell/endothelial/fibroblast neighbourhoods. This project's prespecified `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R` confirms a cohort-wide, quantitative, null-model-calibrated version of the myeloid component specifically: suppressive-myeloid barrier fraction significantly predicts lower clone-tumour engagement (estimate=-0.343, LRT p=0.0069).",
        "novel_contribution": "Cohort-wide statistical confirmation (n=152 clone-sections, calibrated null model, patient-nested mixed model) of a pattern the source paper only illustrates qualitatively/visually.",
    },
    {
        "figure_or_table": "Figure 5",
        "comparison_status": "reproduced_and_extended",
        "this_project_reference": "`06_cell_type_annotation/03_map_external_scrna_reference.py` (map_external_scrna_reference), `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`-12.03 (checkpoint bulk/scRNA projection)",
        "note": "The source paper reports an internal scRNA-to-Xenium label-transfer validation (75.6% RNApred accuracy). This project's External Checkpoint Validation adds an independent external validation (GSE103322, Puram et al. 2017, not the source paper's training data), finding module coherence transfers (4/4 programmes in `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`) even where state-abundance proportions diverge (the Cycling-state discrepancy, `12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`-12.03).",
        "novel_contribution": "An independent (different patient cohort, different sequencing platform) external cross-validation, not merely an internal train/test split on the source paper's data -- a finding that transfers at the gene-module level even where it does not at the abundance level.",
    },
    {
        "figure_or_table": "Figure 6",
        "comparison_status": "methodological_choice",
        "this_project_reference": "Niche and Ecosystem Discovery (niches, composition-vector consensus clustering)",
        "note": "The source paper uses DBSCAN + Moran's I to identify immune aggregates. This project's Niche and Ecosystem Discovery deliberately uses a methodologically distinct approach (composition-vector consensus clustering, Schurch et al. 2020-style) -- a stated methodological choice, not an attempted reproduction of the source paper's specific clustering algorithm.",
        "novel_contribution": "An independent niche-discovery method applied to the same underlying tissue architecture -- if the two methods broadly agree on niche structure this is informal cross-method triangulation, but this project does not claim a formal quantitative comparison of the two clustering algorithms themselves.",
    },
    {
        "figure_or_table": "Figure 7",
        "comparison_status": "reproduced_and_extended",
        "this_project_reference": "Clone Spatial Descriptors (clone-descriptor spatial descriptors), Clone Ecology Confirmatory Models (q2_variance_partition_confirmatory, prespecified)",
        "note": "The source paper illustrates clone-level spatial phenotype segregation in one patient (patient 17, 6 hand-selected clones). This project's Clone Spatial Descriptors/13 extends this to a cohort-wide, rarefaction-normalised, null-model-calibrated, variance-partitioned claim across all 261 clone-sections and 10 patients.",
        "novel_contribution": "Systematic, cohort-wide quantification (identity=20.2%, context=50.4%, patient=29.3% of variance, with bootstrap CIs) of a phenomenon the source paper demonstrates only in a small, hand-selected, illustrative set of clones/patients.",
    },
    {
        "figure_or_table": "Figure 8",
        "comparison_status": "reproduced_and_extended",
        "this_project_reference": "`13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R` (q2_variance_partition_confirmatory, prespecified)",
        "note": "The source paper illustrates patient-level architectural variance using 2 clinically-similar-but-spatially-divergent patients (9 and 12). This project's prespecified `13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R` quantifies the same patient-level variance component across all 10 cohort patients (patient=29.3%, 95% CI 2.7-52.0%).",
        "novel_contribution": "A cohort-wide, quantitative patient-variance-component estimate with uncertainty, directly motivated by (not merely similar to) the source paper's 2-patient illustration.",
    },
    {
        "figure_or_table": "Limitations (Discussion)",
        "comparison_status": "addressed_limitation",
        "this_project_reference": "HPV-Stratified Analysis (the analysis, consolidated/capped), External Validation and Generalisation (this phase -- full external validation)",
        "note": "The source paper's authors explicitly name the small (n=10) spatial cohort as insufficient to fully characterise architectural diversity. This project's HPV-Stratified Analysis (capping HPV contrasts to 1-2 prespecified, with a prospective power simulation showing d>=2.5 required for 80% power) and External Validation and Generalisation (this phase, independent external validation) are a direct methodological response to this same named limitation, not an independently-discovered concern.",
        "novel_contribution": "Methodological infrastructure (the HPV contrast cap, the prospective power simulation, this external-validation phase) built specifically to respond to a limitation the original authors already named -- framed as a response, not a novel discovery of the limitation itself.",
    },
    {
        "figure_or_table": "Neighbourhood signalling (within Discussion/Figure 6 text)",
        "comparison_status": "reproduced_and_extended",
        "this_project_reference": "`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`-14.02 (sender_receiver_pairs, ligand_receptor_database, spatial_scores), `14_spatial_interactions_and_barriers/05_prioritise_testable_interactions.py` (prioritisation)",
        "note": "The source paper's CellChat analysis predicts PVR-TIGIT (tumour-to-T-cell) and CD86-mediated costimulatory (myeloid-to-T-cell) signalling as neighbourhood-specific. This project's independently-designed candidate list (`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`, built without reference to the source paper's CellChat output) includes both pairs (`TIGIT_PVR`, `CD28_CD86`) -- `14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`'s graph-null-model-calibrated spatial test finds both significant, cohort-wide (18/18 and 17/18 sections respectively), and `14_spatial_interactions_and_barriers/05_prioritise_testable_interactions.py`'s prioritisation ranks the CD86/myeloid pair in the top 4 of 14 candidates.",
        "novel_contribution": "Independent rediscovery of both source-paper-named pathways via a different method (spatial graph-null-model calibration vs. CellChat's expression-based inference), plus a cohort-wide effect-size/consistency ranking CellChat's single-method output does not provide.",
    },
    {
        "figure_or_table": "cohort_and_methods.replicate_reproducibility",
        "comparison_status": "diverges",
        "this_project_reference": "`04_quality_control/08_assess_replicate_concordance.R` (04_quality_control/08_assess_replicate_concordance.R)",
        "note": "The source paper reports 1 of 7 replicate pairs discordant in TCR-probe detection (not naming which). This project's `04_quality_control/08_assess_replicate_concordance.R` check, under the same metric (CDR3 probe Jaccard index), finds all 7 pairs at exactly 1.0 -- 0/7 discordant, a direct divergence under a literal same-metric comparison.",
        "novel_contribution": "A methodological insight this project's investigation surfaced: the source paper's stated reason for the discordance (\"differences in tissue content and clonal representation\") suggests a quantitative-representation effect, which a binary Jaccard presence/absence index cannot detect -- explaining, not merely reporting, why the two analyses can legitimately disagree without contradiction. P13 is nonetheless flagged as a convergent soft outlier across three other independent metrics (pseudobulk correlation, transcript-count shift, Moran's I difference), reported as a caveat, not a formal exclusion.",
    },
]


def validate_comparison_rows(rows: list[dict]) -> None:
    """Structural check: every row has all required non-empty fields and
    a valid `comparison_status`."""
    required_fields = [
        "figure_or_table",
        "comparison_status",
        "this_project_reference",
        "note",
        "novel_contribution",
    ]
    for row in rows:
        for field in required_fields:
            if not row.get(field):
                raise PipelineError(
                    f"Comparison row '{row.get('figure_or_table', '<unknown>')}' is missing a non-empty '{field}'."
                )
        if row["comparison_status"] not in VALID_STATUSES:
            raise PipelineError(
                f"Comparison row '{row['figure_or_table']}' has an invalid comparison_status '{row['comparison_status']}'."
            )


def build_source_paper_comparison(project_root: Path) -> dict:
    """Writes Table 8 (source-paper comparison) from `SOURCE_PAPER_COMPARISON`
    above -- already-completed, already-governance-documented comparisons."""
    output_path = project_root / "reports" / "validation" / "source_paper_comparison.tsv"

    validate_comparison_rows(SOURCE_PAPER_COMPARISON)

    result = pd.DataFrame(SOURCE_PAPER_COMPARISON)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    status_counts = result["comparison_status"].value_counts().to_dict()

    return {
        "n_claims_compared": len(result),
        "status_counts": status_counts,
        "output_path": str(output_path),
    }
