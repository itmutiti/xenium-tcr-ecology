# GSE139324 -- Cillo et al. 2020 (Immunity), HNSCC immune-cell single-cell RNA-seq

**Source:** NCBI GEO accession
[GSE139324](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324)
(public since 2019-11-13). Downloaded 2026-07-11 from
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139324/suppl/GSE139324_RAW.tar`.

**Citation:** Cillo AR, Kürten CHL, Tabib T, et al. "Immune Landscape of
Viral- and Carcinogen-Driven Head and Neck Cancer." *Immunity*.
2020;52(1):183-199.e9. PMID:
[31924475](https://pubmed.ncbi.nlm.nih.gov/31924475/).

**Why this dataset (Phase 16.02's real, SECOND independent external
reference, deliberately different from GSE103322/Puram et al. 2017
already used in Phase 12):** re-using GSE103322 for Phase 16's
`cell_state_signature_generalisation` claim (`governance/validation_
plan.tsv`) would not add independent evidence -- it would just
re-derive Phase 12's already-established result on the same
external reference. GSE139324 is an independent HNSCC scRNA-seq cohort:
131,224 cells (CD45+ immune compartment specifically, 10x Genomics
droplet-based, not Smart-seq2) from 63 samples across 26 HPV-/HPV+
HNSCC patients (paired peripheral blood + tumour-infiltrating immune
cells) plus healthy-donor controls -- a substantial, immune-focused
cohort well-suited to this project's T-cell state/checkpoint programme
validation.

**Scope restriction:** only the 26 `HNSCC_<N>_TIL`
(tumour-infiltrating lymphocyte) samples are used for validation -- the
paired `HNSCC_<N>_PBMC` (peripheral blood), `HD_<N>_PBMC` and
`HD_<N>_Tonsil` (healthy-donor control) samples are downloaded (all 63
samples are present in `raw/`) but are a different biological
compartment, not relevant to validating this project's tumour T-cell
state signatures specifically -- left available for potential future
use, not deleted.

**License / reuse:** NCBI GEO public dataset, no access restriction.

**Files:** `GSE139324_RAW.tar` (569,436,160 bytes) extracted to `raw/`
-- 189 files (63 real samples x 3: `_barcodes.tsv.gz`, `_genes.tsv.gz`,
`_matrix.mtx.gz`, standard 10x CellRanger v2 triplet format, Ensembl
gene IDs + symbols). See `checksums.sha256` for the real download's
integrity hash (of the original tar; the extracted `raw/` contents are
reproducible from it).
