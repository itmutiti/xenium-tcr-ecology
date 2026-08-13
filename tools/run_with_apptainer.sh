#!/usr/bin/env bash
# `tools/run_with_apptainer.sh` -- optional, reviewer-friendly Apptainer
# launcher for the xenium-tcr-ecology computational workflow.
#
# This is Mode 1 of two optional Apptainer execution routes (see
# containers/Apptainer.def and docs/execution_manual/EXECUTION_MANUAL.md
# "Reproducibility"): the entire repository runs inside one project-level
# SIF image via explicit bind mounts, with everything else (native
# Conda/Snakemake execution, `snakemake --cores N` directly) completely
# unchanged and unaffected by this script's existence.
#
# Not required to reproduce the analysis. A user without Apptainer can
# still run the complete workflow exactly as documented in README.md /
# EXECUTION_MANUAL.md -- this script is purely additive.
#
# Usage:
#   apptainer build xenium-tcr-ecology.sif containers/Apptainer.def   # once
#   tools/run_with_apptainer.sh snakemake --cores 8                  # any command
#   tools/run_with_apptainer.sh snakemake --cores 8 project_ready     # a specific target
#   tools/run_with_apptainer.sh pytest tests/unit tests/simulation    # anything, not just snakemake
#
# Environment variables (all optional):
#   XENIUM_SIF              Path to the built image (default: xenium-tcr-ecology.sif in the repo root)
#   APPTAINER_CACHEDIR       Apptainer's own build/pull cache (default: Apptainer's own default, ~/.apptainer/cache)
#   APPTAINER_TMPDIR         Apptainer's own temp directory during build (default: $TMPDIR or /tmp)
#   XENIUM_APPTAINER_TMPDIR  Passed into the container as TMPDIR for the pipeline's own temp files (default: /tmp inside the container)
# Override APPTAINER_CACHEDIR/APPTAINER_TMPDIR before building on a host
# with a small $HOME quota (common on shared HPC systems) -- e.g.:
#   export APPTAINER_CACHEDIR=/scratch/$USER/apptainer-cache
set -Eeuo pipefail

info() { echo "[INFO]  $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$HERE")"
[[ -f "$PROJECT_ROOT/manifests/project_paths.yaml" ]] || \
  error "Could not locate the project root from $HERE. Run this script from within the xenium-tcr-ecology checkout."

command -v apptainer >/dev/null 2>&1 || \
  error "apptainer not found on PATH. This is optional -- the native Conda/Snakemake route (README \"Setup\") does not require it."

SIF="${XENIUM_SIF:-$PROJECT_ROOT/xenium-tcr-ecology.sif}"
[[ -f "$SIF" ]] || error "Image not found at '$SIF'. Build it first:
  apptainer build \"$SIF\" \"$PROJECT_ROOT/containers/Apptainer.def\""

[[ $# -gt 0 ]] || error "Usage: $0 <command> [args...]   e.g. $0 snakemake --cores 8"

info "Running inside $SIF (bind: $PROJECT_ROOT -> /workspace, --cleanenv)"
# Wraps every command with `conda run -n xenium-tcr-ecology` explicitly,
# rather than relying on the image's own PATH setup alone -- verified:
# a SIF converted from the Dockerfile-built OCI image (one of
# two supported build routes, see containers/Apptainer.def's header)
# does not expose `snakemake`/`python3`/`Rscript` on a bare `apptainer
# exec`'s PATH (only via its converted Docker ENTRYPOINT, i.e.
# `apptainer run`, or an explicit `conda run` wrapper as used here) --
# this wrapper works correctly regardless of which route built the image.
exec apptainer exec \
  --cleanenv \
  --bind "$PROJECT_ROOT:/workspace" \
  --pwd /workspace \
  --env "TMPDIR=${XENIUM_APPTAINER_TMPDIR:-/tmp}" \
  "$SIF" \
  conda run --no-capture-output -n xenium-tcr-ecology "$@"
