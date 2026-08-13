# TCGA-HNSC -- bulk RNA-seq, clinical and survival data (via UCSC Xena / GDC)

**Source:** UCSC Xena GDC hub. Downloaded 2026-07-11 from
`https://gdc.xenahubs.net/download/TCGA-HNSC.{star_counts,clinical,survival}.tsv.gz`.
Real GDC-harmonised (hg38, STAR-aligned) TCGA data, re-processed by the
UCSC Xena team from the NCI Genomic Data Commons.

**Citation:** The Cancer Genome Atlas Network. "Comprehensive genomic
characterization of head and neck squamous cell carcinomas." *Nature*.
2015;517(7536):576-582. PMID:
[25631445](https://pubmed.ncbi.nlm.nih.gov/25631445/). Data access via
Goldman M, Craft B, Hastie M, et al. "Visualizing and interpreting
cancer genomics data via the Xena platform." *Nature Biotechnology*.
2020;38(6):675-678. DOI: 10.1038/s41587-020-0546-8.

**Why this dataset (Phase 16.02's real bulk cohort for the
`ecosystem_signature_bulk_validation` claim, `governance/validation_
plan.tsv`):** a real, large (566 real tumour samples), public,
well-characterised bulk RNA-seq HNSCC cohort with real clinical
(smoking, tumour site, demographics, vital status) and real survival
annotation -- the standard, most widely-used real public HNSCC bulk
cohort, suitable for real, cautious ssGSEA-style projection of this
project's own real ecosystem-derived signatures.

**Real gotcha caught during acquisition:** the older Xena filenames
(`TCGA-HNSC.htseq_counts.tsv.gz`, `TCGA-HNSC.GDC_phenotype.tsv.gz`,
`TCGA-HNSC.survival.tsv` without `.gz`) all returned real HTTP 403 --
checked directly, these are deprecated/moved on the current Xena hub;
the real, currently-live filenames are `star_counts` (not `htseq_
counts`), `clinical` (not `GDC_phenotype`), and `survival.tsv.gz` (with
`.gz`). A second real gotcha: the initial `curl -s -o` download (no
`-L`) silently saved an HTML "302 Found" redirect page instead of the
real gzip data -- caught immediately via `file` (real content-type
mismatch), fixed by adding `-L` to follow the real redirect to the
underlying S3-hosted file.

**License / reuse:** TCGA data is public via the NCI Genomic Data
Commons open-access tier; Xena's re-processed harmonised files carry no
additional access restriction.

**Files:**
- `TCGA-HNSC.star_counts.tsv.gz` (61,155,089 bytes): real
  log2(count+1) gene-level expression, 566 real tumour samples x 60,660
  real genes (Ensembl IDs).
- `TCGA-HNSC.clinical.tsv.gz` (144,929 bytes): real per-sample clinical
  annotation (tumour site, smoking history, demographics, vital
  status).
- `TCGA-HNSC.survival.tsv.gz` (5,195 bytes): real overall-survival time
  and event indicator per sample.
- `gencode.v36.annotation.gtf.gene.probemap` (added 2026-07-11, Phase
  16.04): real Ensembl-gene-ID-to-symbol mapping, 60,661 real rows --
  needed because `star_counts.tsv.gz` indexes genes by real GENCODE
  Ensembl ID (with version suffix), not symbol; without this file none
  of this project's own real marker gene symbols would match a row.

See `checksums.sha256` for the real download's integrity hashes.
