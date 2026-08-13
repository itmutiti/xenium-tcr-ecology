# references/

Versioned reference tables consumed (not produced-and-consumed-internally)
by phase scripts: marker/cell-type registries (Cell Type Annotation), panel-supported
ligand-receptor interactions (Spatial Interactions and Barriers). Distinct from `governance/`, which
holds this project's own audit trail, and from `data/external/`, which holds
downloaded third-party datasets (External Validation and Generalisation).

Populated by:
- `scripts/06_cell_type_annotation/00_compile_marker_and_reference_registry.py` → `cell_type_marker_registry.tsv`
- `scripts/14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py` → `panel_supported_interactions.tsv`

Empty at scaffold time.
