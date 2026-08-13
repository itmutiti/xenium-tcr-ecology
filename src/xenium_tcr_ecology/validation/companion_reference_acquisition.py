"""Downloads (if not already present) and verifies GSE287301 -- McCord et
al. 2026's own companion scRNA-seq dataset (`data/external/GSE287301/
README.md`): the aggregated gene-expression matrix, used for reference-
label transfer in `06_cell_type_annotation/03_map_external_scrna_
reference.py`, and the 16 per-sample paired scTCR-seq VDJ archives
(`data/external/GSE287301/vdj/README.md`), used as independent ground
truth in `08_tcr_clonal_analysis/09_validate_probe_clones_against_paired_
vdj_ground_truth.py`. Both public, unauthenticated NCBI GEO downloads --
no dedicated acquisition step existed for either before this module; both
consuming scripts previously assumed the files were already staged on
disk with no verification at all.

The GSM-accession-to-chip/pool mapping below is not guessed or inferred
from sequential numbering -- each was confirmed against GEO's
own per-sample record (`Sample_supplementary_file_1` in
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSM>&targ=self&form=text&view=quick`).
"""

from __future__ import annotations

import shutil
import tarfile
import tempfile
from pathlib import Path

from xenium_tcr_ecology.infra.download import download_file, verify_checksums
from xenium_tcr_ecology.infra.exceptions import PipelineError

GEX_SUPPL_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287301/suppl"
GEX_FILES = {
    "GSE287301_filtered_feature_bc_matrix.tar.gz": 2_752_469_184,
    "GSE287301_patient_matrix.txt.gz": 155,
    "GSE287301_aggregation.csv.gz": 570,
}

# GSM accession -> chip/pool directory name, confirmed against
# each sample's own GEO record, not inferred.
VDJ_SAMPLES = {
    "GSM8743474": "chip1pool1",
    "GSM8743475": "chip1pool2",
    "GSM8743476": "chip1pool3",
    "GSM8743477": "chip1pool4",
    "GSM8743478": "chip1pool5",
    "GSM8743479": "chip1pool6",
    "GSM8743480": "chip1pool7",
    "GSM8743481": "chip1pool8",
    "GSM8743482": "chip2pool1",
    "GSM8743483": "chip2pool2",
    "GSM8743484": "chip2pool3",
    "GSM8743485": "chip2pool4",
    "GSM8743486": "chip2pool5",
    "GSM8743487": "chip2pool6",
    "GSM8743488": "chip2pool7",
    "GSM8743489": "chip2pool16",
}
VDJ_MEMBERS = ["clonotypes.csv", "filtered_contig_annotations.csv"]


def _gsm_series_group(gsm: str) -> str:
    """GEO sample suppl directory convention: last 3 digits become 'nnn'."""
    return f"{gsm[:-3]}nnn"


def ensure_gse287301_gex_acquired(project_root: Path) -> Path:
    """Downloads the aggregated gene-expression matrix and its companion
    manifests if not already present at their expected sizes, extracts
    the matrix tarball if `filtered_feature_bc_matrix/` doesn't already
    exist, then verifies every retained file's SHA-256 checksum."""
    dataset_dir = project_root / "data" / "external" / "GSE287301"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    for filename, expected_size in GEX_FILES.items():
        download_file(
            f"{GEX_SUPPL_BASE_URL}/{filename}",
            dataset_dir / filename,
            expected_size_bytes=expected_size,
        )

    matrix_dir = dataset_dir / "filtered_feature_bc_matrix"
    if not matrix_dir.is_dir():
        matrix_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(dataset_dir / "GSE287301_filtered_feature_bc_matrix.tar.gz") as tf:
            tf.extractall(matrix_dir, filter="data")

    checksum_results = verify_checksums(dataset_dir)
    failed = [f for f, ok in checksum_results.items() if not ok]
    if failed:
        raise PipelineError(
            f"Checksum verification failed for {failed} in '{dataset_dir}' -- "
            "re-download required, not silently trusted."
        )
    return dataset_dir


def ensure_gse287301_vdj_acquired(project_root: Path) -> Path:
    """Downloads each of the 16 per-sample VDJ archives to a temporary
    location (not retained -- matching this dataset's existing
    convention, see `data/external/GSE287301/vdj/README.md`) and extracts
    only `clonotypes.csv` and `filtered_contig_annotations.csv` from
    each, skipping any pool whose extracted files already exist. No
    checksums.sha256 is recorded for this subdirectory (by design, see
    that README) -- GEO's own immutable per-sample accessions are the
    provenance chain."""
    vdj_dir = project_root / "data" / "external" / "GSE287301" / "vdj"

    for gsm, pool_name in VDJ_SAMPLES.items():
        pool_dir = vdj_dir / pool_name
        if all((pool_dir / member).is_file() for member in VDJ_MEMBERS):
            continue

        pool_dir.mkdir(parents=True, exist_ok=True)
        url = (
            f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{_gsm_series_group(gsm)}/{gsm}/suppl/"
            f"{gsm}_{pool_name}.tar.gz"
        )
        with tempfile.TemporaryDirectory() as tmp:
            archive_path = Path(tmp) / f"{gsm}_{pool_name}.tar.gz"
            download_file(url, archive_path)
            with tarfile.open(archive_path) as tf:
                for member in VDJ_MEMBERS:
                    tf.extract(f"{pool_name}/{member}", path=Path(tmp), filter="data")
                    # shutil.move, not Path.replace()/os.rename(): the
                    # temporary directory and pool_dir can be on different
                    # filesystems (e.g. Docker's internal /tmp vs. the
                    # /workspace bind mount), which os.rename() cannot
                    # cross (OSError: Invalid cross-device link) -- found
                    # by a real Docker clean-room run, never native
                    # (where both are typically the same filesystem).
                    shutil.move(str(Path(tmp) / pool_name / member), str(pool_dir / member))

    return vdj_dir


def build_companion_reference_acquisition_summary(project_root: Path) -> dict:
    gex_dir = ensure_gse287301_gex_acquired(project_root)
    vdj_dir = ensure_gse287301_vdj_acquired(project_root)

    pools_complete = sum(
        1
        for pool_name in VDJ_SAMPLES.values()
        if all((vdj_dir / pool_name / member).is_file() for member in VDJ_MEMBERS)
    )

    return {
        "gex_dir": str(gex_dir),
        "vdj_dir": str(vdj_dir),
        "vdj_pools_acquired": pools_complete,
        "vdj_pools_expected": len(VDJ_SAMPLES),
    }
