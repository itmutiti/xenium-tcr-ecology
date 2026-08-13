# GSE103322 -- Puram et al. 2017 (Cell), HNSCC single-cell RNA-seq

**Source:** NCBI GEO accession
[GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322)
(public since 2017-11-30). Downloaded 2026-07-11 from
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz`.

**Citation:** Puram SV, Tirosh I, Parikh AS, et al. "Single-Cell
Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in
Head and Neck Cancer." *Cell*. 2017;171(7):1611-1624.e24.
PMID: [29198524](https://pubmed.ncbi.nlm.nih.gov/29198524/).

**Why this dataset (Phase 12's independent external reference):**
5,902 single cells (Smart-seq2 full-length scRNA-seq) from 18 HNSCC
patients across multiple institutions -- independent of this
project's primary data (`GSE287301`, McCord et al. 2026's
companion cohort) and of `GSE287301`'s reference-mapping use in Phase
6.03 (same source-paper cohort, not an independent check). Cell-type
annotations are embedded directly in the supplementary matrix (tumour
site, malignant/non-malignant classification, non-cancer cell type:
Fibroblast, T cell, Endothelial, B cell, Mast, Macrophage, Dendritic,
myocyte) -- 1,237 annotated T cells, sufficient for an
independent check of this project's T-cell-state marker scoring.
Full-transcriptome Smart-seq2 coverage gives complete overlap with every
marker gene this project's cytotoxicity/exhaustion/proliferation/Treg
programs use (`src/xenium_tcr_ecology/preprocess/program_scores.py`,
`src/xenium_tcr_ecology/annotation/t_cell_substates.py`).

**License / reuse:** NCBI GEO public dataset, no access restriction;
original authors' raw FASTQ files are explicitly NOT available ("unable
to provide the raw data due to privacy concerns" per the GEO record) --
only the processed, de-identified expression matrix is used here, which
is exactly what GEO's public supplementary-file mechanism is for.

**File:** `GSE103322_HNSCC_all_data.txt.gz` (90,226,510 bytes). Format: a
plain-text matrix, tab-separated, 23,692 rows -- 1 header row (cell IDs),
5 cell-level metadata rows (`processed by Maxima enzyme`, tumour
site, `classified as cancer cell`, `classified as non-cancer cells`,
`non-cancer cell type`), then 23,686 gene rows (values are the authors'
`log2(TPM/10 + 1)` normalisation, single-quoted gene symbols e.g.
`'GZMB'`). See `checksums.sha256` for the download's integrity hash.
