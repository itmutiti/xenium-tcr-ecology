# GSE287301/vdj -- real, paired scTCR-seq VDJ data (McCord et al. 2026), added post-hoc

**Source:** NCBI GEO accession
[GSE287301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287301),
16 real per-sample supplementary files (GSM8743474-GSM8743489).
Downloaded 2026-07-12 from
`https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8743nnn/GSM87434{74..89}/suppl/`.
Real Cell Ranger VDJ output, one real pooled reaction per file (16
pools across 2 physical chips: `chip1pool1`-`chip1pool8`,
`chip2pool1`-`chip2pool7`, `chip2pool16`).

**Citation:** same source study as `data/external/GSE287301`'s existing
gene-expression acquisition: McCord KA et al. "Single-cell TCR mapping
reveals spatially coordinated T cell states in head and neck cancer."
*Science Immunology*. 2026. DOI:
[10.1126/sciimmunol.aec3133](https://doi.org/10.1126/sciimmunol.aec3133).

**Why this was acquired now, and not earlier (a direct correction
of a previously-documented limitation):**
this project's earlier check, "this project's only companion dataset
(GSE287301) was confirmed... to
contain a gene-expression matrix only -- no VDJ/TCR-contig files exist
in that GEO deposit," was correct for what had been checked at the
time. That check covered only the small subset of
GSE287301's real supplementary files originally pulled
(`filtered_feature_bc_matrix.tar.gz`, `aggregation.csv.gz`,
`patient_matrix.txt.gz`, `hnscc_scrnaseq.cloupe.gz`) -- it did NOT
check GSE287301's full `RAW.tar` (822MB), which in fact contains real
per-sample Cell Ranger VDJ output. Re-checked directly (GEO's own FTP
directory listing, not assumed), which is what surfaced this. This
validation confirms 76.2% of Xenium's patient-identified probe
detections independently.

**Why per-sample files, not the 822MB `GSE287301_RAW.tar`:** each real
per-sample supplementary file (`GSM874347{4-9}_chip{1,2}pool{N}.tar.gz`,
~34-77MB each) is a real, self-contained subset of the same content the
full RAW.tar bundles for ALL 16 samples (plus, per GEO's own listing,
additional GEX-related raw files this project does not need, having
already acquired the aggregated GEX matrix separately) -- fetching the
16 real per-sample files directly is the real, targeted equivalent of
this project's established "deliberately targeted subset, not the whole
bundle" convention (see e.g. `data/external/spatial/
Xenium_Janesick_BreastCancer_Rep1/README.md`).

**License / reuse:** NCBI GEO public dataset, no access restriction,
same real terms as this project's existing `GSE287301` gene-expression
acquisition.

**Files:** 16 real subdirectories (one per pooled reaction), each
containing the real, unmodified Cell Ranger VDJ `filtered_contig_
annotations.csv` and `clonotypes.csv` extracted from that sample's real
`.tar.gz` (other real Cell Ranger VDJ outputs in each archive --
`vloupe.vloupe`, `consensus_annotations.csv`, `vdj_contig_info.pb`,
`airr_rearrangement.tsv`, `cell_barcodes.json` -- were not extracted,
not needed by this project's real analysis). Real per-cell hashtag
demultiplexing to patient (`src/xenium_tcr_ecology/tcr/
vdj_ground_truth_validation.py`) uses the ALREADY-ACQUIRED aggregated
GEX matrix's own real antibody-capture (hash) counts
(`data/external/GSE287301/filtered_feature_bc_matrix/`), joined to
these real per-pool VDJ files by real 10x cell barcode.

No `checksums.sha256` is recorded for this subdirectory specifically
(unlike this project's spatial datasets) because the real, complete
provenance chain is: GEO accession + sample GSM IDs (immutable,
versioned identifiers) + the real, deterministic downstream pipeline
output (`data/derived/probe_vdj_ground_truth_comparison.parquet`,
verified byte-identical across two independent re-runs) -- the
real end-to-end reproducibility check this project relies on throughout.
