"""Clone-level spatial and ecological descriptors (rarefaction-normalised),
discrete-vs-continuum structure testing, variance-partition models, and
barrier-topology / candidate spatial interaction scoring.

The thesis-defining analytical core. Clone Spatial Descriptors
(11.00-11.07) and Clone Ecology Confirmatory Models (13.00-13.06) are
complete. `11_clone_spatial_descriptors/05_test_discrete_vs_continuous_
structure.R`'s prespecified `q2_discrete_vs_continuous_structure_test`
concluded continuous structure -- `taxonomy_version = "v1_provisional"`
(data/releases/v1_clone_structure/) is a 1-factor continuous score, not a
discrete taxonomy. Any future clustering/taxonomy function added here
must not assume discrete categories exist; check the frozen
`structure_type` first.

**Q2 confirmatory result (`13_clone_ecology_confirmatory_models/01_fit_variance_partition_models.R`):** clone
ecological-structure-score variance decomposes into patient 29.3%,
identity 20.2%, context 50.4% -- all three non-trivial, confirming the
prespecified hypothesis. Robust to clone-size adjustment (13.03) and to
leaving out any single patient (13.04, 10/10 folds stable), and to a
segmentation-robustness check (13.05). A standardised atlas of all 261
clone-sections is at `reports/clone_ecology/clone_atlas.html` (13.06).

Spatial Interactions and Barriers is complete. `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s
prespecified `q3_barrier_topology_confirmatory` found a significant
barrier-topology effect on clone-tumour engagement beyond cell-intrinsic
state and niche composition alone.

**Ratio vs. difference normalisation (`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`):** `compute_
engagement_ratio` (tumour_engagement.py) is only valid for bounded
[0, 1] composition-fraction quantities; a centred/z-scored continuous
score (e.g. a gene-programme score) must use `compute_score_excess`
(apc_support.py) instead -- a ratio against a near-zero or negative
baseline is mathematically undefined, not merely imprecise (see
`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`).

Every Clone Spatial Descriptors+ script restricts to `metadata/sample_manifest.tsv`'s
`included_in_primary_hnscc_cohort == True` sections
(`P10_run1`'s specimen is an ameloblastoma, not HNSCC).
"""
