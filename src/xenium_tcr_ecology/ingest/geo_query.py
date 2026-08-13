"""Query GSE300147's GEO FTP filelist.txt and record a snapshot
(`02_raw_data_ingestion/00_query_geo_accession.py`). Public, unauthenticated -- no credentials involved.

GEO series convention (verified against this accession): https://ftp.ncbi.nlm.nih.gov/geo/series/<GSEnnn>/<accession>/suppl/
contains one filelist.txt manifest (Archive/File, Name, Time, Size, Type)
describing everything in that directory, including the single _RAW.tar.

Sends `xenium_tcr_ecology.infra.download.DEFAULT_HEADERS`: NCBI's FTP-
over-HTTPS host returns HTTP 403 for `requests`' default User-Agent and
200 for a descriptive one -- not an authentication issue, but
required for this request to succeed reliably regardless of network path.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import requests

from xenium_tcr_ecology.infra.download import DEFAULT_HEADERS
from xenium_tcr_ecology.infra.exceptions import PipelineError

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def series_group(accession: str) -> str:
    """GEO suppl directory convention: last 3 digits of the accession become 'nnn'."""
    return f"{accession[:-3]}nnn"


def suppl_url(accession: str) -> str:
    return f"{FTP_BASE}/{series_group(accession)}/{accession}/suppl/"


def fetch_filelist(accession: str, timeout: int = 60) -> list[dict]:
    url = suppl_url(accession) + "filelist.txt"
    try:
        resp = requests.get(url, headers=DEFAULT_HEADERS, timeout=timeout)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise PipelineError(f"Failed to fetch '{url}': {exc}") from exc

    lines = resp.text.strip().splitlines()
    if not lines or not lines[0].startswith("#Archive/File"):
        raise PipelineError(
            f"Unexpected filelist.txt format at '{url}': first line was {lines[0]!r}"
        )

    entries = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) != 5:
            raise PipelineError(f"Malformed filelist.txt row in '{url}': {line!r}")
        kind, name, mtime, size, ftype = parts
        entries.append(
            {
                "kind": kind,
                "name": name,
                "remote_mtime": mtime,
                "size_bytes": int(size),
                "file_type": ftype,
            }
        )
    return entries


def build_geo_snapshot(accession: str, output_path: Path) -> dict:
    entries = fetch_filelist(accession)

    archive_entries = [e for e in entries if e["kind"] == "Archive"]
    file_entries = [e for e in entries if e["kind"] == "File"]
    if len(archive_entries) != 1:
        raise PipelineError(
            f"Expected exactly one Archive entry in filelist.txt for {accession}, found {len(archive_entries)}."
        )

    gsm_ids = sorted({e["name"].split("_")[0] for e in file_entries if e["name"].startswith("GSM")})

    snapshot = {
        "accession": accession,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "suppl_url": suppl_url(accession),
        "archive": archive_entries[0],
        "sample_count_from_filelist": len(gsm_ids),
        "gsm_ids": gsm_ids,
        "file_entries": file_entries,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(snapshot, indent=2))

    return {
        "archive_size_bytes": archive_entries[0]["size_bytes"],
        "sample_count": len(gsm_ids),
        "file_count": len(file_entries),
    }
