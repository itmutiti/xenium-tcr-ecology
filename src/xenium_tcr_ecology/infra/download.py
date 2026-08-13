"""Shared HTTP(S) download infrastructure for every external dataset this
project acquires automatically (the primary cohort's GEO archive, and
every external validation/reference dataset under `data/external/`).

**A specific header combination matters, root-caused, not
assumed.** NCBI's GEO FTP-over-HTTPS host (`ftp.ncbi.nlm.nih.gov`) returns
HTTP 403 for Python `requests`' default headers. Isolated by
comparing header-by-header against a working `curl` request (identical
User-Agent, both otherwise default): the differentiator is `requests`'
default `Accept-Encoding: gzip, deflate, br, zstd` -- setting
`Accept-Encoding: identity` (no compression negotiation) alone resolves
it, confirmed both via `requests` with this header and via Python's
stdlib `urllib.request` (which does not send `Accept-Encoding` by
default). Not authentication -- no credential is involved; the same URL,
same file, succeeds with nothing but a different header. Every request
this module makes therefore sends `DEFAULT_HEADERS` by default. 10x
Genomics' CDN and UCSC Xena's hub do not require this, but sending it is
harmless there too.
"""

from __future__ import annotations

import hashlib
import time
from pathlib import Path

import requests

from xenium_tcr_ecology.infra.exceptions import PipelineError

DEFAULT_HEADERS = {
    "User-Agent": (
        "xenium-tcr-ecology-pipeline/1.0 "
        "(https://orcid.org/0009-0006-1768-1887; automated research-data acquisition)"
    ),
    "Accept-Encoding": "identity",
}

_RETRIES = 5
_RETRY_DELAY_SECONDS = 5
_CHUNK_SIZE = 1024 * 1024


def download_file(
    url: str,
    dest: Path,
    *,
    expected_size_bytes: int | None = None,
    headers: dict[str, str] | None = None,
    timeout: int = 60,
) -> None:
    """Streaming download with resume-by-size-check, retry, and a
    descriptive User-Agent. Skips entirely if `dest` already exists and
    matches `expected_size_bytes` (when given) -- safe to call every run,
    re-download only happens when needed."""
    request_headers = {**DEFAULT_HEADERS, **(headers or {})}

    if dest.is_file() and expected_size_bytes is not None:
        if dest.stat().st_size == expected_size_bytes:
            return

    dest.parent.mkdir(parents=True, exist_ok=True)

    last_exc: Exception | None = None
    for attempt in range(1, _RETRIES + 1):
        try:
            with requests.get(
                url, headers=request_headers, stream=True, timeout=timeout, allow_redirects=True
            ) as resp:
                resp.raise_for_status()
                tmp_dest = dest.with_suffix(dest.suffix + ".part")
                with tmp_dest.open("wb") as fh:
                    for chunk in resp.iter_content(chunk_size=_CHUNK_SIZE):
                        if chunk:
                            fh.write(chunk)
                tmp_dest.replace(dest)
            break
        except requests.RequestException as exc:
            last_exc = exc
            if attempt == _RETRIES:
                raise PipelineError(
                    f"Failed to download '{url}' to '{dest}' after {_RETRIES} attempts: {exc}"
                ) from exc
            time.sleep(_RETRY_DELAY_SECONDS)
    else:  # pragma: no cover -- defensive, loop always breaks or raises above
        raise PipelineError(f"Failed to download '{url}': {last_exc}")

    if expected_size_bytes is not None:
        actual = dest.stat().st_size
        if actual != expected_size_bytes:
            dest.unlink(missing_ok=True)
            raise PipelineError(
                f"Downloaded '{url}' to '{dest}' but size ({actual} bytes) does not match "
                f"expected ({expected_size_bytes} bytes) -- transfer incomplete or the "
                "upstream file changed. Re-run to retry."
            )


def verify_checksums(directory: Path) -> dict[str, bool]:
    """Pure, testable: SHA-256 verification of every file listed in
    `directory/checksums.sha256` against its current on-disk content."""
    checksums_path = directory / "checksums.sha256"
    results = {}
    for line in checksums_path.read_text().splitlines():
        if not line.strip():
            continue
        expected_hash, filename = line.split(None, 1)
        filename = filename.strip()
        actual_hash = hashlib.sha256((directory / filename).read_bytes()).hexdigest()
        results[filename] = actual_hash == expected_hash
    return results
