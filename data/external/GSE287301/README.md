# GSE287301 -- McCord et al. 2026's own companion scRNA-seq dataset (gene expression)

**Source:** NCBI GEO accession
[GSE287301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287301).
Downloaded 2026-07-10 from
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287301/suppl/`.

**Citation:** same source study as this project's own primary Xenium
cohort (`data/raw/HNSCC_Xenium/GSE300147`): McCord KA et al. "Single-cell
TCR mapping reveals spatially coordinated T cell states in head and neck
cancer." *Science Immunology*. 2026. DOI:
[10.1126/sciimmunol.aec3133](https://doi.org/10.1126/sciimmunol.aec3133).

**Why this dataset (external reference-label transfer for cell-type
annotation, `06_cell_type_annotation/03_map_external_scrna_reference.py`;
paired scTCR-seq ground truth for TCR probe-clone validation,
`08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_vdj_ground_truth.py`):**
T cells from the same 28-patient HNSCC cohort as this project's own
Xenium data (366,632 cells x 38,606 genes, Cell Ranger 8.0.1) -- found
via a direct search for the source paper's own deposited data.

**License / reuse:** NCBI GEO public dataset, no access restriction.

**Files:**
- `GSE287301_filtered_feature_bc_matrix.tar.gz` (2,752,469,184 bytes):
  aggregated Cell Ranger gene-expression matrix (`matrix.mtx.gz`,
  `barcodes.tsv.gz`, `features.tsv.gz`), extracted to
  `filtered_feature_bc_matrix/`.
- `GSE287301_patient_matrix.txt.gz` (155 bytes): per-sample-to-patient
  mapping.
- `GSE287301_aggregation.csv.gz` (570 bytes): Cell Ranger `aggr`
  aggregation manifest.

See `checksums.sha256` for the real download's integrity hashes.

See `vdj/README.md` for the paired scTCR-seq VDJ data (a separate,
later acquisition of the same GEO series' per-sample supplementary
files).
