"""Archive checksum computation and completeness verification (`02_raw_data_ingestion/02_verify_archive_checksums.py`).

GEO does not publish per-file checksums for this accession (unlike, e.g.,
Dryad) -- there is nothing external to verify the downloaded bytes against.
What this module does instead: (1) compute our own SHA-256 of the archive,
establishing the immutable reference for everything downstream in this
project, and (2) verify the archive's internal file listing (names + sizes,
read via the tarfile module without extracting) exactly matches the
filelist.txt manifest already recorded in metadata/geo_snapshot.json (Phase
2.00) -- catching a truncated or corrupted download before `02_raw_data_ingestion/03_extract_archive_safely.py` spends
time and disk extracting it.
"""

from __future__ import annotations

import hashlib
import json
import tarfile
from pathlib import Path

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter

INTEGRITY_FIELDS = ["filename", "expected_size_bytes", "actual_size_bytes", "status"]


def sha256_of_file(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_archive_and_checksum(
    archive_path: Path,
    geo_snapshot_path: Path,
    sha256sums_path: Path,
    integrity_report_path: Path,
    project_root: Path,
) -> dict:
    if not archive_path.is_file():
        raise PipelineError(
            f"Archive not found: '{archive_path}'. Run `02_raw_data_ingestion/01_download_geo_raw_archive.sh` first."
        )
    if not geo_snapshot_path.is_file():
        raise PipelineError(
            f"Missing '{geo_snapshot_path}'. Run `02_raw_data_ingestion/00_query_geo_accession.py` first."
        )

    snapshot = json.loads(geo_snapshot_path.read_text())
    expected_files = {e["name"]: e["size_bytes"] for e in snapshot["file_entries"]}

    with tarfile.open(archive_path, mode="r") as tf:
        members = {m.name: m.size for m in tf.getmembers() if m.isfile()}

    if integrity_report_path.exists():
        integrity_report_path.unlink()
    writer = InventoryWriter(
        integrity_report_path, project_root=project_root, fields=INTEGRITY_FIELDS
    )

    problems = []
    matched = 0
    for name, expected_size in sorted(expected_files.items()):
        actual_size = members.get(name)
        if actual_size is None:
            status = "missing_from_archive"
            problems.append(f"{name}: missing from archive")
        elif actual_size != expected_size:
            status = "size_mismatch"
            problems.append(f"{name}: expected {expected_size} bytes, archive has {actual_size}")
        else:
            status = "ok"
            matched += 1
        writer.write_row(
            filename=name,
            expected_size_bytes=expected_size,
            actual_size_bytes=actual_size if actual_size is not None else "",
            status=status,
        )

    unexpected = set(members) - set(expected_files)
    for name in sorted(unexpected):
        writer.write_row(
            filename=name,
            expected_size_bytes="",
            actual_size_bytes=members[name],
            status="unexpected_in_archive",
        )
        problems.append(f"{name}: present in archive but not in filelist.txt snapshot")

    archive_hash = sha256_of_file(archive_path)
    sha256sums_path.parent.mkdir(parents=True, exist_ok=True)
    sha256sums_path.write_text(f"{archive_hash}  {archive_path.name}\n")

    if problems:
        raise PipelineError(
            f"Archive integrity check failed ({len(problems)} problem(s)); see '{integrity_report_path}'. "
            f"First few: {problems[:5]}"
        )

    return {
        "archive_sha256": archive_hash,
        "files_verified": matched,
        "files_expected": len(expected_files),
    }
