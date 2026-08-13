"""Barrier topology and candidate spatial interactions.

Backs Spatial Interactions and Barriers (scripts/14_spatial_interactions_and_barriers/). `14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`
(sender-receiver pairs), 14.01 (panel-supported LR pairs), 14.02
(spatially constrained interaction scores), 14.04 (barrier-interface
programme activity) and 14.05 (interaction prioritisation) are
implemented; 14.03 (barrier topology model) is implemented in R
(scripts/14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R, not a
Python module); 14.06 is not yet -- see that stage's scripts before
adding functions here.

Panel gap: `TGF-beta` has no coverage in this project's 623-gene panel
(none of 10 canonical genes are present) -- `sender_receiver_pairs.
TGFB_GENE_SET` is deliberately empty; any script referencing a
"TGF-beta programme" must skip it explicitly rather than approximate it
from unrelated genes.
"""
