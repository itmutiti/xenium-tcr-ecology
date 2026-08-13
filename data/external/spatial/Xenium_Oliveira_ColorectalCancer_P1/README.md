# Xenium_Oliveira_ColorectalCancer_P1 -- de Oliveira et al. 2025 (Nat Genet), 10x Genomics public Xenium colorectal cancer dataset

**Source:** 10x Genomics public dataset (Xenium Onboard Analysis
v2.0.0), sample "P1 Human Colon Cancer". Downloaded 2026-07-12 from
`https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE/`.
Also deposited at NCBI GEO accession
[GSE280314](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280314)
(a Xenium subseries of the study's SuperSeries
[GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318)).

**Citation:** de Oliveira MF, Romero JP, Chung M, Williams SR, Gottscho
AD, Gupta A, Pilipauskas SE, Mohabbat S, Raman N, Sukovich DJ, Patterson
DM, Visium HD Development Team, Taylor SEB. "High-definition spatial
transcriptomic profiling of immune cell populations in colorectal
cancer." *Nature Genetics*. 2025;57(6):1512-1523. DOI:
[10.1038/s41588-025-02193-3](https://doi.org/10.1038/s41588-025-02193-3).

**Why this dataset (second, independent test of the `q1_framework_
generalisation` claim):** a peer-reviewed, single-cell-resolution Xenium
In Situ dataset (the same platform this project's primary cohort and
the existing breast-cancer generalisation test both use) from a
different tissue and cancer type -- colorectal, not
HNSCC, not breast. Chosen over an unpublished 10x preview dataset
(a lung-cancer Immuno-Oncology-panel preview was also identified and
considered) for its peer-reviewed provenance
(Nature Genetics 2025), matching this project's existing standard for
the breast-cancer dataset (Janesick et al., Nature Communications
2023). The Xenium panel used here (10x Human Colon Gene Expression
Panel, 322 genes, supplemented with 100 additional immune-population
genes) was designed by the source study to resolve
diverse immune populations, including T cells -- directly relevant to
this project's T-cell-centric focus, and arguably a better-suited
panel for this generalisation test than the breast-cancer
dataset's general cancer panel.

**Why sample P1 specifically:** the source study profiled 3
patients with Xenium (P1, P2, P5); P1 was chosen as a single
representative section, matching this project's existing
single-representative-section precedent for the breast-cancer dataset
(`Xenium_Janesick_BreastCancer_Rep1`, one of several replicates
in its source study). The `q1_framework_generalisation` claim
concerns the null-model framework's statistical calibration on tissue
topology, not a specific patient's biology, so one
representative section is a sufficient choice, not a statistically
underpowered shortcut.

**License / reuse:** 10x Genomics public datasets are distributed
under CC BY 4.0 (confirmed against 10x Genomics' own public
dataset pages before use).

**Files (a targeted subset of the full ~12.1GB "outs" bundle, not the
whole bundle -- only the analysis outputs this project's validation
actually needs, not the raw microscopy images/transcript-level zarr
archives, mirroring the same targeted-download approach already used
for the breast-cancer dataset):**
- `Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cell_feature_matrix.h5`
  (14,168,777 bytes): per-cell x per-gene count matrix.
- `Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cells.parquet`
  (5,503,158 bytes): per-cell metadata -- `cell_id`, `x_centroid`,
  `y_centroid`, `transcript_counts`, `cell_area`, `nucleus_area`, etc.
  307,762 cells.
- `Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_analysis.tar.gz`
  (51,947,742 bytes, a gzipped tar -- confirmed with `file`
  before extraction, unlike the breast-cancer dataset's mislabelled
  `.tar` which was actually plain POSIX tar despite its filename; this
  file's `.gz` suffix is accurate): extracted to `analysis/`. Graph-based
  unsupervised clustering (`analysis/clustering/
  gene_expression_graphclust/clusters.csv`, 24 clusters, sizes
  504-46,924 cells, 307,383 clustered cells) is used as this
  validation's categorical grouping variable, identical in role to
  the breast-cancer dataset's graph-clustering output.

See `checksums.sha256` for the download's integrity hashes.
