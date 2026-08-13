"""Raw Data Ingestion (02_raw_data_ingestion): Data acquisition, integrity verification and immutable staging."""

rule phase02_00_query_geo_accession:
    """02.00 -- Queries GSE300147 metadata and records accession status, sample list, file sizes, retrieval date and record revision. Primary output: metadata/geo_snapshot.json."""
    input: "results/logs/.checkpoints/project_ready.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/00_query_geo_accession.done")
    shell: "python3 scripts/02_raw_data_ingestion/00_query_geo_accession.py"


rule phase02_01_download_geo_raw_archive:
    """02.01 -- Downloads the raw GEO archive with resumable transfer and conservative retry logic. Primary output: data/raw/GSE300147_RAW.tar."""
    input: "results/logs/02_raw_data_ingestion/.sentinels/00_query_geo_accession.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/01_download_geo_raw_archive.done")
    shell: "bash scripts/02_raw_data_ingestion/01_download_geo_raw_archive.sh"


rule phase02_02_verify_archive_checksums:
    """02.02 -- Calculates SHA-256 checksums, verifies file completeness, produces an immutable manifest. Primary output: data/raw/SHA256SUMS; reports/integrity.tsv."""
    input: "results/logs/02_raw_data_ingestion/.sentinels/01_download_geo_raw_archive.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/02_verify_archive_checksums.done")
    shell: "python3 scripts/02_raw_data_ingestion/02_verify_archive_checksums.py"


rule phase02_03_extract_archive_safely:
    """02.03 -- Extracts without path traversal, preserves timestamps, maps files to GEO sample identifiers. Primary output: data/staged/<sample_id>/."""
    input: "results/logs/02_raw_data_ingestion/.sentinels/02_verify_archive_checksums.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/03_extract_archive_safely.done")
    shell: "python3 scripts/02_raw_data_ingestion/03_extract_archive_safely.py"


rule phase02_04_inventory_xenium_files:
    """02.04 -- Detects H5, Parquet, TIFF and auxiliary files; flags missing mandatory assets and format-version differences. Primary output: metadata/xenium_file_inventory.tsv."""
    input: "results/logs/02_raw_data_ingestion/.sentinels/03_extract_archive_safely.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/04_inventory_xenium_files.done")
    shell: "python3 scripts/02_raw_data_ingestion/04_inventory_xenium_files.py"


rule phase02_05_standardise_sample_directory_layout:
    """02.05 -- Creates a consistent per-section directory structure without modifying source files, using links where possible. Primary output: data/standardised/<section_id>/."""
    input: "results/logs/02_raw_data_ingestion/.sentinels/03_extract_archive_safely.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/05_standardise_sample_directory_layout.done")
    shell: "python3 scripts/02_raw_data_ingestion/05_standardise_sample_directory_layout.py"


rule phase02_06_freeze_raw_data_snapshot:
    """02.06 -- Marks raw and staged data read-only and records provenance and storage usage. Primary output: reports/raw_data_freeze.txt."""
    input:
        "results/logs/02_raw_data_ingestion/.sentinels/04_inventory_xenium_files.done",
        "results/logs/02_raw_data_ingestion/.sentinels/05_standardise_sample_directory_layout.done"
    output: touch("results/logs/02_raw_data_ingestion/.sentinels/06_freeze_raw_data_snapshot.done")
    shell: "bash scripts/02_raw_data_ingestion/06_freeze_raw_data_snapshot.sh"
