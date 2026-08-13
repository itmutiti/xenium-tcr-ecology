"""SpatialData Import (03_spatialdata_import): Xenium import, coordinate harmonisation and canonical data objects."""

rule phase03_00_detect_xenium_format_version:
    """03.00 -- Inspects file schemas, chooses compatible readers, cross-validates against a second independent reader implementation, records deviations. Primary output: metadata/format_versions.tsv."""
    input: "results/logs/.checkpoints/raw_data_frozen.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/00_detect_xenium_format_version.done")
    shell: "python3 scripts/03_spatialdata_import/00_detect_xenium_format_version.py"


rule phase03_01_import_each_section_to_spatialdata:
    """03.01 -- Loads counts, cell metadata, transcript coordinates, boundaries and morphology into one SpatialData object per section. Primary output: data/objects/spatialdata/<section>.zarr."""
    input: "results/logs/.checkpoints/raw_data_frozen.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/01_import_each_section_to_spatialdata.done")
    shell: "python3 scripts/03_spatialdata_import/01_import_each_section_to_spatialdata.py"


rule phase03_02_validate_coordinate_systems:
    """03.02 -- Confirms micrometre units, image orientation, origin and alignment among transcripts, centroids, polygons and images; confirms polygon-vs-centroid convention consistently across sections. Primary output: reports/coordinate_validation/."""
    input:
        "results/logs/03_spatialdata_import/.sentinels/00_detect_xenium_format_version.done",
        "results/logs/03_spatialdata_import/.sentinels/01_import_each_section_to_spatialdata.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/02_validate_coordinate_systems.done")
    shell: "python3 scripts/03_spatialdata_import/02_validate_coordinate_systems.py"


rule phase03_03_create_anndata_expression_objects:
    """03.03 -- Extracts a clean cell-by-gene AnnData object while preserving links to spatial elements. Primary output: data/objects/anndata/<section>.h5ad."""
    input: "results/logs/03_spatialdata_import/.sentinels/01_import_each_section_to_spatialdata.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/03_create_anndata_expression_objects.done")
    shell: "python3 scripts/03_spatialdata_import/03_create_anndata_expression_objects.py"


rule phase03_04_attach_clinical_and_technical_metadata:
    """03.04 -- Joins patient and run metadata using validated keys and prevents silent many-to-many joins. Primary output: updated h5ad/zarr objects; join audit."""
    input: "results/logs/03_spatialdata_import/.sentinels/03_create_anndata_expression_objects.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/04_attach_clinical_and_technical_metadata.done")
    shell: "python3 scripts/03_spatialdata_import/04_attach_clinical_and_technical_metadata.py"


rule phase03_05_build_combined_analysis_object:
    """03.05 -- Concatenates sections with globally unique cell identifiers and explicit patient/run hierarchy. Primary output: data/objects/hnscc_xenium_combined.h5ad."""
    input: "results/logs/03_spatialdata_import/.sentinels/04_attach_clinical_and_technical_metadata.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/05_build_combined_analysis_object.done")
    shell: "python3 scripts/03_spatialdata_import/05_build_combined_analysis_object.py"


rule phase03_06_export_r_interoperability_objects:
    """03.06 -- Creates Seurat-compatible and parquet exports for R analyses without losing spatial coordinates. Primary output: data/objects/r_exports/."""
    input: "results/logs/03_spatialdata_import/.sentinels/04_attach_clinical_and_technical_metadata.done"
    output: touch("results/logs/03_spatialdata_import/.sentinels/06_export_r_interoperability_objects.done")
    shell: "python3 scripts/03_spatialdata_import/06_export_r_interoperability_objects.py"
