# Analysis Amendments

Prespecified analyses and their success criteria are recorded in
`governance/analysis_registry.tsv`. This document records every analysis
added after the initial pipeline was complete: its date, status,
rationale, and relationship to the prespecified analyses. Most entries
below are additive (a sensitivity check, robustness check, external
validation, or presentation change) with no effect on any frozen result.
The most recent entry is the one exception: it corrects which test
statistic gates the release decision for a primary claim, without
changing the claim's success/fail conclusion -- see that entry for detail
and for why it does not fit the "purely additive" pattern of the rest of
this document.

## 2026-07-12 - Cycling-state scoring-method sensitivity check

Script: `12_external_checkpoint_validation/05_rescore_cycling_state_with_primary_method.py`
Status: Completed; hypothesis rejected.
Rationale: The Cycling T-cell state showed the weakest external-reference
concordance of any tested state. This check tested whether the
discrepancy was an artefact of using a simplified scoring proxy for the
external-reference comparison rather than the pipeline's own scoring
method.
Relationship to prespecified analyses: Not prespecified. A post hoc
robustness check using data already available to the pipeline; does not
alter any prespecified claim or its evaluation.

## 2026-07-12 - Q2 variance-partition sensitivity to the Cycling feature

Script: `13_clone_ecology_confirmatory_models/07_test_structure_sensitivity_excluding_cycling.R`
Status: Completed.
Rationale: `cycling_fraction` has the dominant loading on the ecological-
structure score and the weakest external corroboration among the input
features. This check refit the prespecified Q2 variance-partition model
excluding that feature.
Relationship to prespecified analyses: Not prespecified. A post hoc
sensitivity check on the prespecified `q2_variance_partition_confirmatory`
result; found the partition is not robust to excluding this feature.
Reported alongside the prespecified result, not in place of it.

## 2026-07-12 - Q3 barrier-effect covariate-ablation check

Script: `14_spatial_interactions_and_barriers/07_ablate_covariates_for_barrier_effect.R`
Status: Completed.
Rationale: Tested whether the prespecified Q3 barrier-topology
association depends on the full covariate-adjustment set or is driven by
a subset of it.
Relationship to prespecified analyses: Not prespecified. A post hoc
robustness check on the prespecified `q3_barrier_topology_confirmatory`
result; found the association requires the full joint adjustment set.

## 2026-07-12/13 - Second independent spatial dataset (colorectal cancer)

Scripts: `16_external_validation_and_generalisation/08_acquire_second_independent_spatial_dataset.py`,
`16_external_validation_and_generalisation/09_validate_framework_on_second_cancer_type.py`
Status: Completed.
Rationale: Extends the prespecified framework-generalisation claim
(`q1_framework_generalisation`) with a second, independent, peer-reviewed
Xenium dataset in a different tissue type (de Oliveira et al. 2025,
GSE280314).
Relationship to prespecified analyses: Additive external validation of a
prespecified claim, using the same test and success criterion already
applied to the first external dataset. Not itself a new prespecified
analysis; the criterion for `q1_framework_generalisation` was not changed
to accommodate this dataset.

## 2026-07-11 - Freeze-or-revise decision rule for the taxonomy external checkpoint

Script: `12_external_checkpoint_validation/03_decide_freeze_or_revise.py`
Status: Completed.
Rationale: The scaffold specified a predeclared freeze-or-revise rule for
the T-cell-state taxonomy, conditioned on the external module-transfer
and directional-consistency checks (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`, `12_external_checkpoint_validation/02_quantify_directional_consistency.py`). In
this project's execution history, those two checks were implemented in
the same working session as the decision rule itself, before the rule
was finalised in code, so the rule was not literally declared ahead of
seeing the checks' results. Every threshold in the rule is a standard,
generic, field-convention value (`p < 0.05`; a positive sign for rank
correlation), not fitted to reproduce a specific outcome.
Relationship to prespecified analyses: The taxonomy freeze/revise
mechanism itself is not one of the registered primary claims in
`governance/analysis_registry.tsv`; it is pipeline machinery governing
which taxonomy version downstream confirmatory analyses use. Recorded
here for transparency about the rule's timing relative to the checks it
gates.

## 2026-07-13 - Paired scTCR-seq VDJ validation of CDR3-probe patient assignment

Script: `08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_vdj_ground_truth.py`
Status: Completed; updates a previously documented limitation.
Rationale: A documented limitation stated no VDJ/TCR-contig data existed
in the project's companion dataset (GSE287301), based on the file subset
originally acquired. A fuller listing of that dataset's raw archive was
found to contain Cell Ranger VDJ output for all pooled reactions. This
script uses that data to check CDR3-probe-based patient assignment
independently: 76.2% of patient-identified probes were confirmed.
Relationship to prespecified analyses: Not prespecified. Provides
orthogonal, cohort-internal validation of a methodological step
(probe-to-patient assignment) used throughout the TCR pipeline; supersedes
the earlier "no VDJ data available" limitation statement.

## 2026-07-13/21 - Main manuscript figure redesign

Scripts: `17_statistical_closure_and_release/11_build_redesigned_manuscript_figures.py`,
`src/xenium_tcr_ecology/release/redesigned_main_figures.py`
Status: Completed.
Rationale: The original six main figures were finalised before the second
independent spatial dataset and the sensitivity/ablation checks above
existed. The figure set was redesigned to include these findings and to
consolidate two HPV-related figures into one.
Relationship to prespecified analyses: Presentation change, not a
scientific analysis. No result, statistic, or figure source dataset was
recalculated; only which panels are assembled into main-text figures, and
their layout, changed.

## 2026-07-29 - Q3 release-gatekeeping statistic corrected to match the coefficient-level headline claim

Scripts: `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`,
`17_statistical_closure_and_release/01_control_multiplicity_and_report_effects.R`
Status: Completed. Corrects which test statistic gates the release
decision for `q3_barrier_topology_confirmatory`; does not change the
PASS conclusion.

Rationale: `governance/analysis_registry.tsv`'s registered hypothesis for
Q3 is written jointly ("Fibroblast/myeloid barrier topology... explains
residual variance"), and `03_model_barrier_topology_by_structure.R`
already computes the correct statistic for that joint hypothesis -- a
2-degree-of-freedom likelihood-ratio test (LRT) comparing state+niche
against state+niche+fibroblast+myeloid (chisq=9.961, df=2, p=0.0069,
unchanged by this entry, still reported only in this script's own log and
figure, not persisted to disk). However, the paper's headline result and
`barrier_topology_model_results.parquet`'s frozen effect estimate/CI are
about `suppressive_myeloid_barrier_fraction` specifically, a narrower
claim the joint 2-df LRT does not test. Release gatekeeping
(`01_control_multiplicity_and_report_effects.R`) had been recomputing an
ad hoc two-sided Wald p-value from that coefficient's frozen
estimate/se (raw p=0.00152, Bonferroni x5=0.00762) to stand in for a
single-coefficient test that had never actually been computed. This
correction adds a properly nested 1-df LRT specifically for
`suppressive_myeloid_barrier_fraction` (state+niche+fibroblast vs.
+myeloid; chisq=9.082, df=1, raw p=0.00258, Bonferroni x5=0.01291),
persisted as new `lrt_chisq`/`lrt_df`/`lrt_pvalue` columns in
`barrier_topology_model_results.parquet` (populated for the myeloid row
only; existing `estimate`/`se`/`ci_low`/`ci_high` columns are unchanged).
`01_control_multiplicity_and_report_effects.R` now reads and
Bonferroni-corrects this persisted LRT p-value directly, rather than
recomputing a Wald substitute. All three statistics -- the previous Wald
substitute, the new myeloid-only LRT, and the unchanged joint LRT -- pass
Bonferroni correction at alpha=0.05 with comfortable margin; the
release's Q3 PASS conclusion is unchanged.

Because `barrier_topology_model_results.parquet` is one of the files
`00_freeze_primary_results.py` freezes into `data/releases/final_primary/`,
this correction required a deliberate re-freeze: the file's hash changed
(pre-fix `barrier_topology_model_results.parquet` sha256:
`fa8d838bc6fc247f401f73f9825330b3f8550f9f82e76919c0de8f63f70926ac`),
correctly tripping `build_primary_results_freeze()`'s
"refusing to silently re-freeze over a changed scientific result" guard.
The prior `MANIFEST.json`/`checksums.sha256` were replaced (not merged)
by the re-freeze; this note is the record of what changed and why, since
the replaced manifest no longer carries that history itself.

Relationship to prespecified analyses: Does not change
`q3_barrier_topology_confirmatory`'s registered hypothesis, unit of
analysis, or success criterion, and does not change its PASS outcome.
Changes which already-computed-elsewhere statistic the release's
multiplicity-correction step reads for this claim, from a substitute
computed ad hoc in that step to a purpose-built statistic computed and
persisted in the claim's own source script. The joint 2-df LRT that
directly answers the registered hypothesis as literally written is
unchanged and was not affected by this correction.

---

Prespecified analyses and their frozen success criteria are unaffected by
every entry above except the last, which is explicitly a correction to
which statistic gates a release decision (not to the criterion or its
outcome) -- see that entry for detail. See
`governance/analysis_registry.tsv` for full detail on each.
