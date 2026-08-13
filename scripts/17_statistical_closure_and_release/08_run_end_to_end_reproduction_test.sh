#!/usr/bin/env bash
# `17_statistical_closure_and_release/08_run_end_to_end_reproduction_test.sh` -- 08_run_end_to_end_reproduction_test.sh
#
# Smoke test, not a full clean-room reproduction. Runs a representative
# subset of fast, non-external-data-dependent scripts end-to-end
# (matching the scaffold's "first on a test subset" stage) -- the "then
# the full dataset" stage is explicitly not attempted here, a deliberate
# scope decision (this project's full pipeline took many hours across
# development; a full-dataset re-run is not a reasonable undertaking to
# verify within a single session). See
# src/xenium_tcr_ecology/release/reproduction_test.py's module
# docstring for
# the reproducibility gap (missing r-arrow/r-testthat in
# environment.lock) this exact process found and fixed.
#
# Primary output: reports/reproduction_test/

set -Eeuo pipefail

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
[[ -n "$PROJECT_ROOT" ]] || { echo "[ERROR] Could not locate project root." >&2; exit 1; }

cd "$PROJECT_ROOT"
python3 -c "
from pathlib import Path
from xenium_tcr_ecology.release.reproduction_test import build_reproduction_test

summary = build_reproduction_test(Path('.'))
print(f\"[OK]   {summary['n_pass']}/{summary['n_scripts_tested']} smoke-test script(s) PASS.\")
print(f\"[OK]   {summary['scope_note']}\")
print(f\"[OK]   Wrote {summary['output_path']}\")
if summary['n_fail'] > 0:
    raise SystemExit(1)
"

# Meta-check on the DAG itself, not a DAG step: every rule in
# workflow/rules/*.smk declares only a sentinel touch-file as its
# Snakemake `output:` -- the actual scientific artefact is named only in
# the rule's docstring, so Snakemake's dependency graph cannot notice if
# a script's output silently failed to appear. This closes that gap for
# whatever sentinels exist at this point in the run.
python3 tools/verify_declared_outputs.py --project-root "$PROJECT_ROOT"
