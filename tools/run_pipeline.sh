#!/usr/bin/env bash
# `tools/run_pipeline.sh` -- reviewer-facing launcher for the complete
# xenium-tcr-ecology computational workflow. Selects, in strict priority
# order, the first usable execution route:
#
#   1. Docker     (preferred)
#   2. Apptainer  (fallback)
#   3. Native Conda/Mamba + Snakemake (final fallback -- always supported,
#      see README.md "Setup")
#
# All three routes invoke the identical Snakefile, config, and rules, and
# all three build from the identical, committed
# environment/conda/environment.lock (the canonical, byte-pinned
# dependency lock -- NOT gitignored; see that file's own directory and
# EXECUTION_MANUAL.md "Environment setup" for why committing it is safe
# for the containers specifically even though the native route still
# resolves environment/conda/main.yml fresh on each arbitrary host).
# Backend selection never depends on environment.lock, any other
# generated file, or any machine-specific record -- Docker's own
# usability is checked first, purely as a capability probe (daemon
# reachable, a container actually runs), with zero reference to what
# has or hasn't been generated locally. This script selects and prepares
# the software environment only, never the science. See
# docs/execution_manual/EXECUTION_MANUAL.md "Reproducibility" for what
# has and has not been verified for each route, and CHANGELOG.md / git
# history for how each fix here was found.
#
# Usage:
#   tools/run_pipeline.sh                              # auto backend, snakemake --cores all
#   tools/run_pipeline.sh snakemake --cores 8 project_ready
#   tools/run_pipeline.sh --backend docker snakemake --cores 8
#   tools/run_pipeline.sh --backend apptainer pytest tests/unit tests/simulation
#
# --backend auto|docker|apptainer|native (default: auto, or
# $XENIUM_FORCE_ROUTE if set):
#   - auto: probes Docker, then Apptainer, then native, in that order,
#     and runs the pipeline through the first one that is usable. Every
#     rejected backend's reason is printed and logged;
#     which backend was ultimately selected is always printed and
#     logged clearly.
#   - docker | apptainer | native: uses that backend specifically and
#     NEVER silently falls back to another one. If the requested backend
#     is not usable, this script explains why, states explicitly that no
#     automatic fallback has been performed because a specific backend
#     was requested, detects which alternative backends are available
#     on this host, and prints the exact command to use each
#     one (or to return to automatic selection) before exiting non-zero.
#
# Fallback (in --backend auto only) is conservative and transparent: a
# route is skipped only for a verified environment/capability failure
# discovered *before* its command starts running -- no Docker daemon,
# insufficient Docker permissions, a build that fails, unsupported
# Apptainer host configuration, or failure to create the native
# environment. Once the selected route's command has actually started,
# this script does not fall back further if it fails -- it aborts and
# reports that route's own error, so a scientific run never silently
# continues on a different route with partially-generated, mixed-
# provenance outputs.
#
# Selected backend, every rejected backend's reason, the container image
# tag and digest (or Apptainer .sif path and identifier), and the
# dependency-lock checksum are always printed and saved to
# results/logs/run_pipeline/route_selection.log and
# results/logs/run_pipeline/provenance.json.
#
# Environment variables (all optional):
#   XENIUM_DOCKER_IMAGE         Local Docker tag to use/build
#                                (default: xenium-tcr-ecology:local)
#   XENIUM_DOCKER_IMAGE_DIGEST  Immutable published image reference, e.g.
#                                ghcr.io/<owner>/<repo>@sha256:<digest> --
#                                if set, pulled and preferred over any
#                                local build. Not yet set by this
#                                repository -- see the release-time
#                                finalisation checklist in
#                                docs/execution_manual/EXECUTION_MANUAL.md.
#   XENIUM_SIF                  Path to a pre-built .sif
#                                (default: xenium-tcr-ecology.sif in the repo root)
#   XENIUM_APPTAINER_TMPDIR     TMPDIR exported inside the Apptainer container
#   XENIUM_FORCE_ROUTE          Legacy synonym for --backend (docker |
#                                apptainer | native); an explicit --backend
#                                flag takes precedence over this.
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$HERE")"
[[ -f "$PROJECT_ROOT/manifests/project_paths.yaml" ]] || {
    echo "[run_pipeline] ERROR: could not locate the project root from $HERE. Run this script from within the xenium-tcr-ecology checkout." >&2
    exit 1
}
cd "$PROJECT_ROOT"

LOG_DIR="$PROJECT_ROOT/results/logs/run_pipeline"
mkdir -p "$LOG_DIR"
ROUTE_LOG="$LOG_DIR/route_selection.log"
PROVENANCE_JSON="$LOG_DIR/provenance.json"
: >"$ROUTE_LOG"

info() { echo "[run_pipeline] $*" | tee -a "$ROUTE_LOG"; }
error() {
    echo "[run_pipeline] ERROR: $*" | tee -a "$ROUTE_LOG" >&2
    exit 1
}

lock_checksum() {
    if [[ -f "$PROJECT_ROOT/environment/conda/environment.lock" ]]; then
        sha256sum "$PROJECT_ROOT/environment/conda/environment.lock" | cut -d' ' -f1
    else
        echo "absent"
    fi
}

record_route() {
    local backend="$1" identifier="$2" digest="$3"
    echo "$backend" >"$LOG_DIR/selected_route"
    info "=== Backend selected: $backend ($identifier) ==="
    cat >"$PROVENANCE_JSON" <<JSON
{
  "backend": "$backend",
  "identifier": "$identifier",
  "digest": "$digest",
  "environment_lock_sha256": "$(lock_checksum)",
  "command": "${CMD[*]}",
  "timestamp_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
    info "Provenance written to $PROVENANCE_JSON"
}

# ------------------------------------------------------------ Arg parsing --
BACKEND="${XENIUM_FORCE_ROUTE:-auto}"
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
    --backend)
        BACKEND="$2"
        shift 2
        ;;
    --backend=*)
        BACKEND="${1#--backend=}"
        shift
        ;;
    *)
        ARGS+=("$1")
        shift
        ;;
    esac
done
case "$BACKEND" in
auto | docker | apptainer | native) ;;
*) error "Unknown backend '$BACKEND' (from --backend or \$XENIUM_FORCE_ROUTE; expected auto, docker, apptainer, or native)" ;;
esac

if [[ ${#ARGS[@]} -eq 0 ]]; then
    CMD=(snakemake --cores all)
else
    CMD=("${ARGS[@]}")
fi

# ---------------------------------------------------------------- Probing --
# Capability checks only -- no build, no execution. Used both to decide
# auto-mode fallback and to report available alternatives when
# an explicitly-requested backend is rejected.
DOCKER_REJECT_REASON=""
probe_docker() {
    command -v docker >/dev/null 2>&1 || {
        DOCKER_REJECT_REASON="docker client not found on PATH"
        return 1
    }
    docker info >/dev/null 2>&1 || {
        DOCKER_REJECT_REASON="daemon not reachable (not running, or this user lacks permission -- e.g. not in the 'docker' group)"
        return 1
    }
    docker run --rm hello-world >/dev/null 2>&1 || {
        DOCKER_REJECT_REASON="a minimal container run failed"
        return 1
    }
    return 0
}

APPTAINER_REJECT_REASON=""
probe_apptainer() {
    command -v apptainer >/dev/null 2>&1 || {
        APPTAINER_REJECT_REASON="apptainer not found on PATH"
        return 1
    }
    return 0
}

NATIVE_REJECT_REASON=""
probe_native() {
    command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1 || {
        NATIVE_REJECT_REASON="neither mamba nor conda found on PATH"
        return 1
    }
    return 0
}

print_alternatives() {
    info "Detecting which alternative backends are available on this host..."
    local any=0
    if probe_docker; then
        info "  - docker (available):    tools/run_pipeline.sh --backend docker ${CMD[*]}"
        any=1
    fi
    if probe_apptainer; then
        info "  - apptainer (available): tools/run_pipeline.sh --backend apptainer ${CMD[*]}"
        any=1
    fi
    if probe_native; then
        info "  - native (available):    tools/run_pipeline.sh --backend native ${CMD[*]}"
        any=1
    fi
    [[ "$any" -eq 1 ]] || info "  (no execution backend is currently usable on this host -- see the reasons logged above)"
    info "Or return to automatic backend selection: tools/run_pipeline.sh ${CMD[*]}"
}

reject_explicit_backend() {
    info "Backend '$BACKEND' was explicitly requested but is not usable: $1"
    info "No automatic fallback has been performed because a specific backend was requested."
    print_alternatives
    exit 1
}

# ------------------------------------------------------------------ Docker --
run_docker() {
    local image="${XENIUM_DOCKER_IMAGE:-xenium-tcr-ecology:local}"
    if [[ -n "${XENIUM_DOCKER_IMAGE_DIGEST:-}" ]]; then
        info "Docker: pulling verified image $XENIUM_DOCKER_IMAGE_DIGEST"
        if docker pull "$XENIUM_DOCKER_IMAGE_DIGEST" >>"$ROUTE_LOG" 2>&1; then
            image="$XENIUM_DOCKER_IMAGE_DIGEST"
        else
            info "Docker: pull of verified image failed"
            return 1
        fi
    elif docker image inspect "$image" >/dev/null 2>&1; then
        info "Docker: using existing local image '$image'"
    else
        # Deliberately no check for environment/conda/environment.lock
        # here: it is normally committed and present, but even if it
        # were absent, containers/Dockerfile bootstraps a
        # fresh resolve from environment/conda/main.yml automatically --
        # backend selection must never depend on this file, and it
        # doesn't.
        info "Docker: no local image found; building (environment preparation, not scientific execution)"
        if ! docker build -t "$image" -f containers/Dockerfile "$PROJECT_ROOT" >>"$ROUTE_LOG" 2>&1; then
            info "Docker: image build failed (see $ROUTE_LOG)"
            return 1
        fi
    fi

    local digest
    digest="$(docker image inspect "$image" --format '{{.Id}}' 2>/dev/null || echo unknown)"
    record_route "docker" "$image" "$digest"
    exec docker run --rm \
        --user "$(id -u):$(id -g)" \
        -v "$PROJECT_ROOT:/workspace" \
        "$image" \
        "${CMD[@]}"
}

# --------------------------------------------------------------- Apptainer --
run_apptainer() {
    local sif="${XENIUM_SIF:-$PROJECT_ROOT/xenium-tcr-ecology.sif}"
    if [[ -f "$sif" ]]; then
        info "Apptainer: using existing image '$sif'"
    else
        info "Apptainer: no .sif found; building from containers/Apptainer.def (environment preparation, not scientific execution)"
        if apptainer build --fakeroot "$sif" "$PROJECT_ROOT/containers/Apptainer.def" >>"$ROUTE_LOG" 2>&1; then
            : # built directly
        elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
            info "Apptainer: --fakeroot build unavailable on this host (see EXECUTION_MANUAL.md 'Reproducibility'); trying docker-daemon:// conversion"
            if ! {
                docker build -t xenium-tcr-ecology:for-sif -f containers/Dockerfile "$PROJECT_ROOT" >>"$ROUTE_LOG" 2>&1 \
                    && apptainer build "$sif" docker-daemon://xenium-tcr-ecology:for-sif:latest >>"$ROUTE_LOG" 2>&1
            }; then
                info "Apptainer: docker-daemon:// conversion also failed (see $ROUTE_LOG)"
                return 1
            fi
        else
            info "Apptainer: cannot build (--fakeroot unavailable on this host and no Docker available to convert from)"
            return 1
        fi
    fi

    local sif_id
    sif_id="$(sha256sum "$sif" 2>/dev/null | cut -d' ' -f1 || echo unknown)"
    record_route "apptainer" "$sif" "$sif_id"
    exec apptainer exec \
        --cleanenv \
        --bind "$PROJECT_ROOT:/workspace" \
        --pwd /workspace \
        --env "TMPDIR=${XENIUM_APPTAINER_TMPDIR:-/tmp}" \
        "$sif" \
        conda run --no-capture-output -n xenium-tcr-ecology "${CMD[@]}"
}

# ------------------------------------------------------------------ Native --
run_native() {
    local mgr
    mgr="$(command -v mamba || command -v conda)"

    if ! "$mgr" env list 2>/dev/null | grep -qw "xenium-tcr-ecology"; then
        info "Native: environment 'xenium-tcr-ecology' not found; creating from environment/conda/main.yml"
        "$mgr" env create -f "$PROJECT_ROOT/environment/conda/main.yml" >>"$ROUTE_LOG" 2>&1 \
            || error "Native: environment creation failed -- see $ROUTE_LOG"
    fi

    local smk
    smk="$(command -v snakemake || true)"
    if [[ -z "$smk" && -x "$PROJECT_ROOT/.venv/bin/snakemake" ]]; then
        smk="$PROJECT_ROOT/.venv/bin/snakemake"
    fi
    if [[ -z "$smk" ]]; then
        info "Native: snakemake not found on PATH or in .venv/; creating the documented orchestration .venv"
        python3.13 -m venv "$PROJECT_ROOT/.venv" >>"$ROUTE_LOG" 2>&1 \
            && "$PROJECT_ROOT/.venv/bin/pip" install snakemake >>"$ROUTE_LOG" 2>&1 \
            || error "Native: .venv/snakemake setup failed -- see $ROUTE_LOG"
        smk="$PROJECT_ROOT/.venv/bin/snakemake"
    fi

    record_route "native" "$mgr env: xenium-tcr-ecology" "n/a"
    # Only a bare `snakemake` invocation is substituted with the resolved
    # binary -- other commands (pytest, Rscript, ...) run exactly as
    # given, as a reviewer who activated the conda env by hand would run
    # them.
    if [[ "${CMD[0]}" == "snakemake" ]]; then
        CMD=("$smk" "${CMD[@]:1}")
    fi
    # No --no-capture-output here: `mamba run` (unlike `conda run`, used
    # unconditionally inside the Docker/Apptainer routes above, where the
    # binary is always genuine conda) breaks on this flag -- confirmed
    # ("exec: --: invalid option" from mamba's own generated
    # wrapper script) -- and $mgr may resolve to either.
    exec "$mgr" run -n xenium-tcr-ecology "${CMD[@]}"
}

# ----------------------------------------------------------------- Select --
case "$BACKEND" in
docker)
    info "Backend: docker (explicit -- no automatic fallback)"
    if probe_docker; then
        run_docker
        error "Docker was requested and passed capability checks, but failed during build/execution (see $ROUTE_LOG). No automatic fallback has been performed because a specific backend was requested."
    fi
    reject_explicit_backend "$DOCKER_REJECT_REASON"
    ;;
apptainer)
    info "Backend: apptainer (explicit -- no automatic fallback)"
    if probe_apptainer; then
        run_apptainer
        error "Apptainer was requested and passed capability checks, but failed during build/execution (see $ROUTE_LOG). No automatic fallback has been performed because a specific backend was requested."
    fi
    reject_explicit_backend "$APPTAINER_REJECT_REASON"
    ;;
native)
    info "Backend: native (explicit -- no automatic fallback)"
    if probe_native; then
        run_native
    fi
    reject_explicit_backend "$NATIVE_REJECT_REASON"
    ;;
auto)
    info "Backend: auto -- probing Docker -> Apptainer -> native in priority order"
    info "Command: ${CMD[*]}"
    if probe_docker; then
        run_docker || info "Docker: passed capability checks but build/execution failed -- falling back"
    else
        info "Docker: rejected -- $DOCKER_REJECT_REASON"
    fi
    if probe_apptainer; then
        run_apptainer || info "Apptainer: passed capability checks but build/execution failed -- falling back"
    else
        info "Apptainer: rejected -- $APPTAINER_REJECT_REASON"
    fi
    if probe_native; then
        run_native
    else
        error "Native: rejected -- $NATIVE_REJECT_REASON. This is the final fallback; no backend is usable on this host."
    fi
    ;;
esac
