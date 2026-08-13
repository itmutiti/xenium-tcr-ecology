# Xenium_Janesick_BreastCancer_Rep1 -- Janesick et al. 2023 (Nat Commun), 10x Genomics public Xenium breast cancer dataset

**Source:** 10x Genomics public dataset (production pipeline v1.0.1),
"In Situ Sample 1, Replicate 1". Downloaded 2026-07-11 from
`https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/`.
Also deposited at NCBI GEO accession
[GSE243280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243280)
(the SuperSeries cited in the paper's own Data Availability statement),
specifically its Xenium SubSeries
[GSE243168](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243168),
sample GSM7780153 ("Breast Cancer, Xenium In Situ Spatial Gene
Expression Rep 1"). An earlier version of this README, and this
project's dataset registry, recorded GSE243275 -- this is GEO's own
scRNA-seq/Visium SubSeries of the same SuperSeries and contains no
Xenium samples; corrected 2026-07-19 after direct verification against
the live GEO record.

**Citation:** Janesick A, Shelansky R, Gottscho AD, et al. "High
resolution mapping of the tumor microenvironment using integrated
single-cell, spatial and in situ analysis." *Nature Communications*.
2023;14:8353. DOI:
[10.1038/s41467-023-43458-x](https://doi.org/10.1038/s41467-023-43458-x).

**Why this dataset (Phase 16.01's independent test of Q1 --
framework generalisation, `governance/validation_plan.tsv`'s
`q1_framework_generalisation` claim):** single-cell resolution
Xenium In Situ data (the same platform/technology this project's
primary cohort uses) from an independent cohort (FFPE
human breast cancer, not HNSCC, not McCord et al.'s GSE300147
cohort) -- the field's standard, most widely-used public Xenium
benchmark dataset (extensively reused across the spatial transcriptomics
literature, e.g. cell-type-annotation benchmarking studies). The
Q1 claim this dataset tests is about the framework's statistical
calibration (Type I error, power of the null-model calibration suite,
Phase 9.08) -- a platform/methods-level question, not a claim
about HNSCC biology specifically, so a different cancer type
on the same single-cell-resolution platform is the correct
choice (an HNSCC-specific spatial dataset with per-cell
clonal/lineage information at Xenium's single-cell resolution does
not appear to exist publicly at time of acquisition; considered
alternatives included GSE252265, a Visium/spot-level HNSCC dataset, but
not single-cell resolution).

**License / reuse:** 10x Genomics public datasets are distributed under
CC BY 4.0 (confirmed against 10x Genomics' own public dataset pages
before use).

**Files (a targeted subset of the full ~9.86GB "outs" bundle, not the
whole bundle -- only the analysis outputs this project's validation
actually needs, not the raw microscopy images/transcript-level zarr
archives):**
- `Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5`
  (12,148,885 bytes): per-cell x per-gene count matrix.
- `Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet` (3,453,894
  bytes): per-cell metadata -- `cell_id`, `x_centroid`,
  `y_centroid`, `transcript_counts`, `cell_area`, `nucleus_area`, etc.
  167,780 cells.
- `Xenium_FFPE_Human_Breast_Cancer_Rep1_analysis.tar` (64,317,440
  bytes, a plain POSIX tar despite the upstream `.gz`-suffixed filename
  -- confirmed with `file` before extraction failed once and was
  corrected): extracted to `analysis/`. Graph-based unsupervised
  clustering
  (`analysis/clustering/gene_expression_graphclust/clusters.csv`, 20
  clusters, sizes 400-24,022 cells) is used as this validation's
  categorical "lineage" grouping variable -- Q1's claim
  concerns framework calibration (whether the null-model machinery is
  correctly calibrated for any categorical spatial labelling), so
  an unsupervised cluster label is a valid, sufficient grouping
  variable; biologically naming each cluster (e.g. "epithelial",
  "fibroblast") is not required for this claim and was not attempted
  here.

See `checksums.sha256` for the download's integrity hashes.
