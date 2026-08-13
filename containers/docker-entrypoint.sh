#!/usr/bin/env bash
# Entrypoint for xenium-tcr-ecology's Docker image.
#
# Reviewers are expected to run this image with `--user "$(id -u):$(id -g)"`
# (documented in EXECUTION_MANUAL.md and used by tools/run_pipeline.sh) so
# that files written to the bind-mounted repository are owned by the host
# user, not root. An arbitrary host UID has no /etc/passwd entry inside the
# container by default, which breaks anything that calls getpwuid()
# (Python's getpass.getuser(), used indirectly by some dependencies) or
# needs a writable $HOME (numba's function cache, matplotlib, R). This
# script self-registers the running UID/GID before handing off, the same
# pattern used by Red Hat's and Jupyter's "arbitrary UID"-compatible
# images. /etc/passwd and /etc/group are made group-writable at image
# build time (see Dockerfile) specifically to allow this.
set -Eeuo pipefail

if ! getent passwd "$(id -u)" >/dev/null 2>&1; then
    echo "xenium-runtime-user:x:$(id -u):$(id -g)::/tmp:/bin/bash" >> /etc/passwd
fi
if ! getent group "$(id -g)" >/dev/null 2>&1; then
    echo "xenium-runtime-group:x:$(id -g):" >> /etc/group
fi

export HOME=/tmp

# Single-threaded BLAS by default (overridable), matching
# containers/Apptainer.def's %environment -- Docker was missing this
# entirely. Multi-threaded BLAS makes floating-point summation order
# depend on core count and thread scheduling, which is not just a
# performance question: it measurably changes results that depend on
# iterative numerical refitting under an otherwise-fixed RNG seed (found
# directly -- a bootstrap-based confirmatory model estimate differed by
# ~7% between a lower-core-count host and a 16-core Docker container,
# despite an identical seed, because the seed governs the random draws,
# not the arithmetic each refit performs on them).
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"

exec conda run --no-capture-output -n xenium-tcr-ecology "$@"
