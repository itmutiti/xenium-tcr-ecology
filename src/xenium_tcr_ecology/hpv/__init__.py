"""HPV-stratified analysis (consolidated, capped)
(scripts/15_hpv_stratified_analysis/).

`config/reproducibility_policy.yaml` (`gates.hpv_single_contrast`) and
`governance/analysis_registry.tsv`'s `hpv_primary_contrast` reserved slot
cap this project to 1-2 prespecified confirmatory HPV contrasts, fixed by
`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`
before any test is run -- a deliberate constraint against small-sample
(n~5 vs n~5) HPV subgroup fishing.

All 7 milestones are implemented. `15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py` (HPV status
validation), 15.01 (the single primary HPV contrast, n=4 vs n=4, fixed
in `governance/hpv_primary_contrasts.yaml`) and 15.06 (the final claim-
strength synthesis, `results/hpv_claim_strength.tsv`) are Python;
15.02-15.05 (power simulation, composition/structure comparisons,
robustness checks) are implemented in R under `scripts/15_hpv_stratified_analysis/`.
"""
