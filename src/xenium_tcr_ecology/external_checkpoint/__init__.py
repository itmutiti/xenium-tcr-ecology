"""Early external sanity-check and taxonomy-freeze gate.

Backs External Checkpoint Validation (scripts/12_external_checkpoint_validation/). All 5 milestones
(12.00-12.03, 12.05) are implemented. 12.03 applies a predeclared, fully
automated decision rule (freeze iff all programs transfer AND overall
Spearman rho > 0; revise otherwise) and writes the result directly to
`governance/freeze_decision.tsv` -- no separate attribution or approval
step exists in this DAG.

**Outcome: freeze.** `taxonomy_version` remains `v1_provisional`
(governance/freeze_decision.tsv) -- the discrepancy between this
project's T-cell-state proportions and the independent reference is
concentrated in 2 of 7 states (`Cycling` -- `11_clone_spatial_descriptors/06_discover_provisional_structure.R`'s single
dominant continuous-structure-axis driver -- and `Ambiguous`), flagged as
a carried-forward caveat; the other 5 states show 80% pairwise
directional agreement, and all 4 transcriptional programmes (including
the one driving `Cycling`) transfer externally as coherent gene modules
.

**Indexing note (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`):** `t_cell_states.parquet`'s index is a plain
positional `RangeIndex`, not cell IDs -- always `.set_index("cell_id")`
(or select the `cell_id` column directly) before joining/matching
against it. An earlier version of `program_transfer.py` used `.index`
directly and matched zero cells against `primary_analysis_matrix.h5ad`,
running its computation on an empty dataframe without erroring.
"""
