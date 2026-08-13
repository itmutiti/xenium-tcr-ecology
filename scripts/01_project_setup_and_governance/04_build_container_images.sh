#!/usr/bin/env bash
# `01_project_setup_and_governance/04_build_container_images.sh` -- 04_build_container_images.sh
#
# Builds the Docker image from containers/Dockerfile, which itself builds
# from the byte-pinned environment/conda/environment.lock (`01_project_setup_and_governance/03_lock_software_environments.sh`'s
# output) -- not a re-resolve of main.yml. Refuses to run if the lock file
# doesn't exist yet, rather than silently building from an unpinned spec.
#
# Apptainer/.sif conversion (for HPC use) is not attempted here: it
# requires either a running Apptainer install or a remote build service,
# neither of which is assumed present. Convert the built Docker image with
# `apptainer build release.sif docker-daemon://xenium-tcr-ecology:latest`
# when Apptainer is available and HPC deployment is actually needed --
# fabricating a .sif build step now, unvalidated, would be worse than
# leaving it as a documented manual step.
#
# Primary output: a local Docker image tagged xenium-tcr-ecology:latest
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

command -v docker >/dev/null 2>&1 || error "docker not found on PATH."
docker info >/dev/null 2>&1 || error "docker daemon not reachable (is it running?)."

LOCK_FILE="${PROJECT_ROOT}/environment/conda/environment.lock"
[[ -f "$LOCK_FILE" ]] || error "Missing ${LOCK_FILE}. Run scripts/01_project_setup_and_governance/03_lock_software_environments.sh first."

IMAGE_TAG="xenium-tcr-ecology:latest"
info "Building ${IMAGE_TAG} from containers/Dockerfile (context: ${PROJECT_ROOT})..."
docker build -f "${PROJECT_ROOT}/containers/Dockerfile" -t "$IMAGE_TAG" "$PROJECT_ROOT"

ok "Built ${IMAGE_TAG}."
info "Smoke-testing the image (import the package + spatial stack inside it)..."
docker run --rm "$IMAGE_TAG" python -c "import xenium_tcr_ecology, scanpy, squidpy, spatialdata; print('container smoke test OK, xenium_tcr_ecology', xenium_tcr_ecology.__version__)"

ok "Image smoke test passed."
info "Apptainer .sif conversion not run automatically -- see this script's header comment. Run manually if/when HPC deployment is needed:"
info "  apptainer build ${PROJECT_ROOT}/containers/release.sif docker-daemon://${IMAGE_TAG}"
