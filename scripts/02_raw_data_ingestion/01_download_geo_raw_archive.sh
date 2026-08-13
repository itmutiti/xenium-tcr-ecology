#!/usr/bin/env bash
# `02_raw_data_ingestion/01_download_geo_raw_archive.sh` -- 01_download_geo_raw_archive.sh
#
# Downloads GSE300147_RAW.tar (~52 GiB) from the GEO FTP suppl directory,
# public/unauthenticated. Resumable via curl's native -C - (byte-range
# resume) and --retry with backoff, so a transfer failure retries from
# where it left off rather than restarting from zero -- a known failure
# mode on large single files.
#
# Verifies the final file size against the archive size recorded in
# metadata/geo_snapshot.json (`02_raw_data_ingestion/00_query_geo_accession.py`) -- refuses to declare success on
# a truncated file even if curl itself exits 0.
#
# Sends the same headers as xenium_tcr_ecology.infra.download.DEFAULT_HEADERS
# (a descriptive User-Agent, Accept-Encoding: identity): NCBI's host 403s
# Python `requests`' default headers (root-caused directly to the
# Accept-Encoding negotiation, not authentication -- see that module's
# docstring). curl's own defaults do not trigger this, but the same
# headers are sent here too for consistency and future robustness.
#
# Primary output: data/raw/GSE300147_RAW.tar
set -Eeuo pipefail

ACCESSION="GSE300147"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/${ACCESSION}/suppl/${ACCESSION}_RAW.tar"
USER_AGENT="xenium-tcr-ecology-pipeline/1.0 (https://orcid.org/0009-0006-1768-1887; automated research-data acquisition)"

info() { echo "[INFO]  $*"; }
ok() { echo "[OK]    $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

PROJECT_ROOT="${1:-}"
if [[ -z "$PROJECT_ROOT" && -n "${XENIUM_TCR_ECOLOGY_ROOT:-}" ]]; then
  PROJECT_ROOT="$XENIUM_TCR_ECOLOGY_ROOT"
fi
if [[ -z "$PROJECT_ROOT" ]]; then
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  candidate="$here"
  while [[ "$candidate" != "/" ]]; do
    if [[ -f "$candidate/manifests/project_paths.yaml" ]]; then
      PROJECT_ROOT="$candidate"
      break
    fi
    candidate="$(dirname "$candidate")"
  done
fi
[[ -n "$PROJECT_ROOT" ]] || error "Could not locate project root. Pass it as \$1 or set \$XENIUM_TCR_ECOLOGY_ROOT."

command -v curl >/dev/null 2>&1 || error "curl not found on PATH."
command -v python3 >/dev/null 2>&1 || error "python3 not found on PATH."

SNAPSHOT="${PROJECT_ROOT}/metadata/geo_snapshot.json"
[[ -f "$SNAPSHOT" ]] || error "Missing ${SNAPSHOT}. Run scripts/02_raw_data_ingestion/00_query_geo_accession.py first."

EXPECTED_SIZE="$(python3 -c "import json; print(json.load(open('${SNAPSHOT}'))['archive']['size_bytes'])")"
[[ -n "$EXPECTED_SIZE" ]] || error "Could not read expected archive size from ${SNAPSHOT}."

DEST_DIR="${PROJECT_ROOT}/data/raw"
DEST="${DEST_DIR}/${ACCESSION}_RAW.tar"
mkdir -p "$DEST_DIR"

if [[ -f "$DEST" ]]; then
  CURRENT_SIZE="$(stat -c%s "$DEST" 2>/dev/null || stat -f%z "$DEST")"
  if [[ "$CURRENT_SIZE" -eq "$EXPECTED_SIZE" ]]; then
    ok "${DEST} already present and matches expected size (${EXPECTED_SIZE} bytes). Skipping download."
    exit 0
  fi
  info "Existing partial file ($CURRENT_SIZE / $EXPECTED_SIZE bytes) -- resuming."
fi

info "Downloading ${URL}"
info "Expected size: ${EXPECTED_SIZE} bytes ($(python3 -c "print(f'{${EXPECTED_SIZE}/1e9:.2f} GB')"))"

curl -C - --retry 10 --retry-delay 5 --retry-all-errors --fail --show-error \
  -A "$USER_AGENT" -H "Accept-Encoding: identity" \
  -o "$DEST" "$URL"

FINAL_SIZE="$(stat -c%s "$DEST" 2>/dev/null || stat -f%z "$DEST")"
if [[ "$FINAL_SIZE" -ne "$EXPECTED_SIZE" ]]; then
  error "Downloaded file size ($FINAL_SIZE bytes) does not match expected ($EXPECTED_SIZE bytes) -- transfer incomplete or archive changed on the server. Re-run this script to resume."
fi

ok "Downloaded ${DEST} (${FINAL_SIZE} bytes, matches expected size)."
