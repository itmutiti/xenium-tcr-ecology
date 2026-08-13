# HNSCC Xenium Spatial TCR Ecology

[![CI](https://github.com/itmutiti/xenium-tcr-ecology/actions/workflows/ci.yml/badge.svg)](https://github.com/itmutiti/xenium-tcr-ecology/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A calibrated, generalisable statistical and software framework for
clone-resolved spatial ecology, developed on the public GSE300147 Xenium
HNSCC cohort (McCord et al., *Science Immunology* 2026) and validated
against five further independent external datasets (see
`manifests/dataset_registry.yaml`).

If you use this software, please cite it - see `CITATION.cff`.

## Installation

```bash
tools/run_pipeline.sh
```

Tries Docker, then Apptainer, then native conda/Snakemake, and runs the
complete pipeline unattended on any Linux host from a bare clone - no
manual environment setup. Use `--backend docker|apptainer|native` to
pick a route explicitly instead of auto-probing.

To install natively instead of through the launcher:

```bash
mamba env create -f environment/conda/main.yml
conda activate xenium-tcr-ecology
pip install -e ".[dev]"
pytest tests/unit tests/simulation
```

Snakemake orchestration runs from its own `.venv/` (Python 3.13), kept
separate from the conda env above:

```bash
python3.13 -m venv .venv
.venv/bin/pip install snakemake
conda activate xenium-tcr-ecology
.venv/bin/snakemake --cores 8 -n     # dry run
.venv/bin/snakemake --cores 8        # full computational analysis
```

If Python 3.13 isn't already available, install it via conda-forge
(`mamba create -n py313_bootstrap python=3.13`) rather than the system
package manager.

## Usage

```bash
snakemake --cores 8            # full computational analysis (default target)
snakemake --cores 8 <target>   # a specific named checkpoint
snakemake --cores 8 publication_release   # public data bundle + code/DOI archive
```

All datasets are downloaded and checksum-verified automatically; no
manual step is required. Expect on the order of hours and ~100 GB of
disk. Outputs land in `figures/`, `tables/`, `reports/`, `results/`,
`data/`, and `manuscript/manuscript.pdf`.

## Documentation

Full documentation - execution order, datasets, computational
requirements, reproducibility verification, architecture, and phase
reference - is at
[docs/execution_manual/EXECUTION_MANUAL.md](docs/execution_manual/EXECUTION_MANUAL.md)
and the [mkdocs site](docs/index.md).

## License and citation

MIT - see `LICENSE`. See `CITATION.cff` for how to cite this software;
version, release date, and repository URL are filled in at the first
tagged release.
