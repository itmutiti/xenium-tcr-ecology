"""Downloads (if not already present) and verifies the two external
references Phase 16.02 acquires
(`16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`):
`GSE139324` (Cillo et al. 2020, Immunity -- a second independent HNSCC
scRNA-seq reference, deliberately different from GSE103322/Puram et al.
2017 already used in External Checkpoint Validation) and `TCGA-HNSC`
(bulk RNA-seq + clinical/survival, via UCSC Xena). Both are public,
unauthenticated downloads -- see `data/external/scrna/GSE139324/README.md`
and `data/external/bulk/TCGA-HNSC/README.md` for the full acquisition
provenance, citations and license (the same READMEs this module's exact
URLs and file sizes are taken from).
"""

from __future__ import annotations

import re
import tarfile
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.download import download_file, verify_checksums
from xenium_tcr_ecology.infra.exceptions import PipelineError

TIL_SAMPLE_PATTERN = re.compile(r"(HNSCC_\d+_TIL)")

GSE139324_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139324/suppl/GSE139324_RAW.tar"
)
GSE139324_EXPECTED_SIZE_BYTES = 569_436_160

TCGA_XENA_BASE = "https://gdc.xenahubs.net/download"
TCGA_FILES = {
    "TCGA-HNSC.star_counts.tsv.gz": 61_155_089,
    "TCGA-HNSC.clinical.tsv.gz": 144_929,
    "TCGA-HNSC.survival.tsv.gz": 5_195,
    "gencode.v36.annotation.gtf.gene.probemap": 3_206_945,
}


def find_til_sample_names(filenames: list[str]) -> set[str]:
    """Pure, testable: unique `HNSCC_<N>_TIL` sample names found in a
    list of raw filenames -- excludes PBMC/HD/Tonsil files, a
    deliberate scope restriction to the tumour-infiltrating compartment
    (see module docstring / README)."""
    names = set()
    for filename in filenames:
        match = TIL_SAMPLE_PATTERN.search(filename)
        if match:
            names.add(match.group(1))
    return names


def ensure_gse139324_acquired(project_root: Path) -> Path:
    """Downloads `GSE139324_RAW.tar` if not already present at the
    expected size, extracts it to `raw/` if that directory doesn't
    already exist, then verifies the tar's SHA-256 checksum. Safe to
    call every run."""
    gse139324_dir = project_root / "data" / "external" / "scrna" / "GSE139324"
    tar_path = gse139324_dir / "GSE139324_RAW.tar"
    raw_dir = gse139324_dir / "raw"

    download_file(GSE139324_URL, tar_path, expected_size_bytes=GSE139324_EXPECTED_SIZE_BYTES)

    if not raw_dir.is_dir():
        raw_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path) as tf:
            tf.extractall(raw_dir, filter="data")

    return gse139324_dir


def ensure_tcga_hnsc_acquired(project_root: Path) -> Path:
    """Downloads every TCGA-HNSC file (UCSC Xena GDC hub, public,
    unauthenticated -- redirects to S3, `download_file`'s
    `allow_redirects=True` follows this, matching the `-L` gotcha
    documented in `data/external/bulk/TCGA-HNSC/README.md`) if not
    already present at its expected size. Safe to call every run."""
    tcga_dir = project_root / "data" / "external" / "bulk" / "TCGA-HNSC"
    for filename, expected_size in TCGA_FILES.items():
        download_file(
            f"{TCGA_XENA_BASE}/{filename}", tcga_dir / filename, expected_size_bytes=expected_size
        )
    return tcga_dir


def build_scrna_reference_acquisition_summary(project_root: Path) -> dict:
    gse139324_dir = ensure_gse139324_acquired(project_root)
    tcga_dir = ensure_tcga_hnsc_acquired(project_root)
    raw_dir = gse139324_dir / "raw"

    gse139324_checksums = verify_checksums(gse139324_dir)
    tcga_checksums = verify_checksums(tcga_dir)
    failed = [f for f, ok in {**gse139324_checksums, **tcga_checksums}.items() if not ok]
    if failed:
        raise PipelineError(
            f"Checksum verification failed for {failed} -- re-download required, not silently trusted."
        )

    til_samples = find_til_sample_names([p.name for p in raw_dir.iterdir()])
    if len(til_samples) == 0:
        raise PipelineError(f"No HNSCC_<N>_TIL sample files found in '{raw_dir}'.")

    counts = pd.read_csv(tcga_dir / "TCGA-HNSC.star_counts.tsv.gz", sep="\t", nrows=0)
    n_tcga_samples = len(counts.columns) - 1  # first column is Ensembl_ID

    return {
        "gse139324_n_til_samples": len(til_samples),
        "gse139324_n_files_verified": len(gse139324_checksums),
        "tcga_n_samples": n_tcga_samples,
        "tcga_n_files_verified": len(tcga_checksums),
    }
