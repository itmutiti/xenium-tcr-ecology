"""Safe tar extraction, grouped by GEO sample identifier (`02_raw_data_ingestion/03_extract_archive_safely.py`).

Explicit path-traversal validation on every member before extraction
(absolute paths and ".." components rejected), on top of tarfile's own
'data' extraction filter where available (Python >= 3.11.4/3.12) --
defense in depth rather than relying on either alone. Preserves original
mtimes. Extraction is a pure structural operation: it does not decompress
the individual .gz members inside (that happens, if ever, when a later
phase actually reads a file), so staged size stays close to archive size.
"""

from __future__ import annotations

import os
import re
import tarfile
from pathlib import Path

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter

EXTRACTION_FIELDS = ["gsm_accession", "filename", "destination_path", "size_bytes", "status"]

GSM_PREFIX_RE = re.compile(r"^(GSM\d+)_")


def _validate_member_path(name: str) -> None:
    if name.startswith("/") or name.startswith("\\"):
        raise PipelineError(f"Refusing to extract member with absolute path: '{name}'")
    if ".." in Path(name).parts:
        raise PipelineError(f"Refusing to extract member with path-traversal component: '{name}'")


def _gsm_for(name: str) -> str:
    match = GSM_PREFIX_RE.match(name)
    if not match:
        raise PipelineError(f"Could not determine GSM sample id from filename: '{name}'")
    return match.group(1)


def safe_extract(archive_path: Path, dest_root: Path, project_root: Path) -> dict:
    if not archive_path.is_file():
        raise PipelineError(
            f"Archive not found: '{archive_path}'. Run `02_raw_data_ingestion/01_download_geo_raw_archive.sh` first."
        )

    dest_root.mkdir(parents=True, exist_ok=True)
    inventory_path = (
        project_root / "results" / "tables" / "02_raw_data_ingestion" / "extraction_inventory.tsv"
    )
    if inventory_path.exists():
        inventory_path.unlink()
    writer = InventoryWriter(inventory_path, project_root=project_root, fields=EXTRACTION_FIELDS)

    extracted = 0
    with tarfile.open(archive_path, mode="r") as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        for member in members:
            # Defence in depth: reject suspicious member names outright
            # (raises) even though the destination path below is always
            # built from the basename only, never from member.name
            # directly, so a ".." component could not otherwise escape
            # dest_root -- failing loudly on an unexpected name is better
            # than silently stripping it and extracting anyway.
            _validate_member_path(member.name)
            gsm = _gsm_for(Path(member.name).name)
            dest_dir = dest_root / gsm
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest_path = dest_dir / Path(member.name).name

            fileobj = tf.extractfile(member)
            if fileobj is None:
                raise PipelineError(f"Could not read member '{member.name}' from archive.")
            with dest_path.open("wb") as out:
                while chunk := fileobj.read(1 << 20):
                    out.write(chunk)

            os.utime(dest_path, (member.mtime, member.mtime))

            writer.write_row(
                gsm_accession=gsm,
                filename=Path(member.name).name,
                destination_path=dest_path,
                size_bytes=member.size,
                status="extracted",
            )
            extracted += 1

    # dest_root can contain non-sample entries (e.g. the scaffold's own
    # .gitkeep placeholder for an otherwise-empty data/ directory) -- count
    # only sample subdirectories, not every directory entry.
    sample_count = len([p for p in dest_root.iterdir() if p.is_dir()])
    return {"files_extracted": extracted, "samples": sample_count}
