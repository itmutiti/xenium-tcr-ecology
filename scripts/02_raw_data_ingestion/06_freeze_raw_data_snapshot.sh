#!/usr/bin/env bash
# `02_raw_data_ingestion/06_freeze_raw_data_snapshot.sh` -- 06_freeze_raw_data_snapshot.sh
#
# Marks data/raw/ and data/staged/ read-only (chmod a-w, recursive) and
# records a provenance report: archive SHA-256 (from `02_raw_data_ingestion/02_verify_archive_checksums.py`), git
# commit, freeze timestamp, and storage usage per data layer. This is the
# last Raw Data Ingestion script and the phase's own completion gate ("raw data are
# read-only and backed up") is what this enforces.
#
# Reversible (chmod u+w undoes the read-only marking) but deliberately not
# run automatically as part of any other script -- freezing is a discrete,
# visible action taken once `02_raw_data_ingestion/00_query_geo_accession.py`-2.05 have all succeeded, not a side
# effect.
#
# Primary output: reports/raw_data_freeze.txt
set -Eeuo pipefail

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

RAW_DIR="${PROJECT_ROOT}/data/raw"
STAGED_DIR="${PROJECT_ROOT}/data/staged"
SHA256SUMS="${RAW_DIR}/SHA256SUMS"
XENIUM_INVENTORY="${PROJECT_ROOT}/metadata/xenium_file_inventory.tsv"
REPORT="${PROJECT_ROOT}/reports/raw_data_freeze.txt"

[[ -d "$RAW_DIR" ]] || error "Missing ${RAW_DIR}. Run `02_raw_data_ingestion/01_download_geo_raw_archive.sh` first."
[[ -f "$SHA256SUMS" ]] || error "Missing ${SHA256SUMS}. Run `02_raw_data_ingestion/02_verify_archive_checksums.py` first."
[[ -d "$STAGED_DIR" ]] || error "Missing ${STAGED_DIR}. Run `02_raw_data_ingestion/03_extract_archive_safely.py` first."
[[ -f "$XENIUM_INVENTORY" ]] || error "Missing ${XENIUM_INVENTORY}. Run `02_raw_data_ingestion/04_inventory_xenium_files.py` first."

info "Marking ${RAW_DIR} and ${STAGED_DIR} read-only..."
chmod -R a-w "$RAW_DIR"
chmod -R a-w "$STAGED_DIR"

mkdir -p "$(dirname "$REPORT")"
{
  echo "# Raw data freeze report"
  echo "Frozen at (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if GIT_COMMIT="$(git -C "$PROJECT_ROOT" rev-parse HEAD 2>/dev/null)"; then
    echo "Git commit: ${GIT_COMMIT}"
  else
    echo "Git commit: none yet (repository initialised but no commits made)"
  fi
  echo ""
  echo "## Archive checksum"
  cat "$SHA256SUMS"
  echo ""
  echo "## Storage usage"
  echo "data/raw:          $(du -sh "$RAW_DIR" 2>/dev/null | cut -f1)"
  echo "data/staged:       $(du -sh "$STAGED_DIR" 2>/dev/null | cut -f1)"
  echo ""
  echo "## Permissions after freeze"
  echo "data/raw:    $(stat -c '%A' "$RAW_DIR")"
  echo "data/staged: $(stat -c '%A' "$STAGED_DIR")"
} > "$REPORT"

ok "Froze ${RAW_DIR} and ${STAGED_DIR} read-only. Wrote ${REPORT}"
info "To reverse (e.g. to re-run an earlier phase): chmod -R u+w ${RAW_DIR} ${STAGED_DIR}"
