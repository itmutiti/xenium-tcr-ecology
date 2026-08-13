#!/usr/bin/env bash
# `01_project_setup_and_governance/03_lock_software_environments.sh` -- 03_lock_software_environments.sh
#
# Exports the already-built xenium-tcr-ecology conda environment
# (Python + R) to machine-generated lock files, and captures R sessionInfo()
# and OS/compiler details. Does not build the environment itself -- that is
# `mamba env create -f environment/conda/main.yml`, run once when the
# environment doesn't yet exist; this script only locks what already exists,
# matching the human-edited-spec-vs-machine-generated-lock split documented
# in environment/conda/main.yml's header.
#
# Primary output: environment/conda/environment.lock;
#                  environment/conda/environment.resolved.yml;
#                  environment/R_sessionInfo.txt
set -Eeuo pipefail

ENV_NAME="xenium-tcr-ecology"

info() { echo "[INFO]  $*"; }
ok() { echo "[OK]    $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

# Root resolution mirrors xenium_tcr_ecology.infra.paths.find_project_root:
# --project-root -> $XENIUM_TCR_ECOLOGY_ROOT -> walk up to manifests/project_paths.yaml.
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

command -v conda >/dev/null 2>&1 || error "conda not found on PATH."
conda env list | grep -qE "^${ENV_NAME}[[:space:]]" || error "Conda environment '${ENV_NAME}' does not exist yet. Run: mamba env create -f ${PROJECT_ROOT}/environment/conda/main.yml"

ENV_DIR="${PROJECT_ROOT}/environment"
mkdir -p "${ENV_DIR}/conda"

info "Exporting full resolved environment (with build strings) -> environment/conda/environment.resolved.yml"
conda env export -n "$ENV_NAME" > "${ENV_DIR}/conda/environment.resolved.yml"

info "Exporting explicit package+version pin list -> environment/conda/environment.lock"
conda list -n "$ENV_NAME" --explicit --md5 > "${ENV_DIR}/conda/environment.lock" 2>/dev/null \
  || conda list -n "$ENV_NAME" > "${ENV_DIR}/conda/environment.lock"

info "Capturing R sessionInfo() -> environment/R_sessionInfo.txt"
{
  echo "# Captured $(date -u +%Y-%m-%dT%H:%M:%SZ) from conda environment '${ENV_NAME}'"
  echo "# OS: $(uname -a)"
  echo "# gcc: $(conda run -n "$ENV_NAME" gcc --version 2>/dev/null | head -1 || echo 'not found in env')"
  echo ""
  conda run -n "$ENV_NAME" Rscript -e 'sessionInfo()'
} > "${ENV_DIR}/R_sessionInfo.txt"

PY_VERSION="$(conda run -n "$ENV_NAME" python --version 2>&1)"
R_VERSION="$(conda run -n "$ENV_NAME" Rscript -e 'cat(R.version.string)' 2>&1)"

ok "Locked. ${PY_VERSION}; ${R_VERSION}."
ok "Wrote environment/conda/environment.resolved.yml, environment/conda/environment.lock, environment/R_sessionInfo.txt"
