# xenium-tcr-ecology

A calibrated, generalisable statistical and software framework for
clone-resolved spatial ecology, developed on the public GSE300147 Xenium
HNSCC cohort and validated beyond it.

Start with:

- [execution_manual/EXECUTION_MANUAL.md](execution_manual/EXECUTION_MANUAL.md) - execution order,
  datasets, and reproduction instructions.
- [Architecture](architecture.md) - why the repository is shaped this way.
- [Phase reference](phases.md) - the computational workflow stage reference.

**Status:** complete. All 17 computational workflow stages are
implemented and tested. The default Snakemake target covers the
computational analysis only; a separate `publication_release` target
builds the public data bundle and code/DOI archive.
