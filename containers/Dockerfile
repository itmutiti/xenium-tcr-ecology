# Frozen-release container image, built from the byte-pinned, committed
# environment/conda/environment.lock (a `conda list --explicit --md5`
# export -- platform-pinned to linux-64, the only platform this image
# ever builds on, so committing it does not carry the cross-platform-
# drift risk the native route's own fresh `main.yml` resolve has on an
# arbitrary reviewer host; see EXECUTION_MANUAL.md "Environment setup").
# Built by scripts/01_project_setup_and_governance/04_build_container_images.sh
# and re-verified by scripts/17_statistical_closure_and_release/09_run_
# null_model_calibration_regression.py before every tagged release.
#
# Backend selection (tools/run_pipeline.sh) must never depend on this
# file's presence -- environment.lock is now committed and normally
# always present, but if it is ever genuinely absent (e.g. a shallow
# checkout, or deliberately deleted to regenerate it), the build below
# falls back to a fresh `mamba env create -f main.yml` resolve
# automatically rather than failing outright, so Docker still comes up
# in one command; that fallback build is not byte-pinned the way the
# committed lock is, and is logged as such.
#
# Two real gaps found and fixed here (architecture assessment for
# containers/Apptainer.def, which mirrors both fixes -- see that file):
# (1) `environment.lock` is a `conda list --explicit` export, which only
# ever captures conda-managed packages -- `markdown`, `xhtml2pdf`, and
# `build` are declared in pyproject.toml (not conda-installed in the
# real, currently-working dev environment despite also being listed,
# redundantly, in environment/conda/main.yml), so a previous `pip
# install --no-deps -e .` here silently omitted them, breaking
# manuscript rendering and software packaging inside the image. Fixed by
# installing the full documented `.[dev]` extra, matching the native
# route exactly (README "Setup"), not a container-specific shortcut.
# (2) `pdftocairo` (poppler-utils) is required for PNG/SVG figure
# export but was never actually present in the conda environment on the
# development machine -- it worked there only because the host OS
# happened to have system poppler-utils already installed via apt,
# never guaranteed by this project's own setup. Installed explicitly as
# a system package below instead of relying on conda's `poppler`
# package (declared in main.yml but not reliably resolved).
#
# Pinned to an immutable digest, not :latest -- containers/Apptainer.def
# pins the identical digest, so both routes provably start from the same
# base OS layer. Found by direct A/B testing (see docs/reproducibility_
# investigation.md): a fresh Docker build and a fresh Apptainer build,
# same day, same machine, same committed environment.lock, produced a
# ~7% different bootstrap-based confirmatory model estimate; both
# built `FROM ...:latest`, which is not reproducible across builds taken
# at different times if the upstream tag moves in between. Regenerate
# this digest deliberately (`docker pull condaforge/miniforge3:latest`,
# copy the new digest here and into Apptainer.def identically) rather
# than reverting to a floating tag.
FROM condaforge/miniforge3@sha256:8edc74dee5b29511a9c1a8941c77efbf438b47edffb6320ab8675dc06bd618c2

RUN apt-get update && apt-get install -y --no-install-recommends poppler-utils \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# The `environment.loc[k]` glob is deliberate: it matches
# environment.lock when present (the normal, committed case) without
# making COPY fail when it is absent -- main.yml always matches, so the
# COPY itself never errors either way.
COPY environment/conda/main.yml environment/conda/environment.loc[k] /tmp/
RUN if [ -f /tmp/environment.lock ]; then \
        echo "[build] Using committed environment/conda/environment.lock (canonical, byte-pinned)."; \
        conda create --name xenium-tcr-ecology --file /tmp/environment.lock; \
    else \
        echo "[build] environment.lock not found in build context -- bootstrapping from environment/conda/main.yml (fresh resolve, NOT byte-pinned; this should not happen on a normal checkout, where the lock is committed)."; \
        mamba env create --name xenium-tcr-ecology -f /tmp/main.yml; \
    fi && conda clean -afy

COPY pyproject.toml README.md /workspace/
COPY src /workspace/src

SHELL ["conda", "run", "-n", "xenium-tcr-ecology", "/bin/bash", "-c"]
RUN pip install -e ".[dev]"

COPY . /workspace

# Reviewers run this image with `--user "$(id -u):$(id -g)"` (see
# EXECUTION_MANUAL.md and tools/run_pipeline.sh) so that files written to
# the bind-mounted repository are owned by the host user, not root. An
# arbitrary host UID has no /etc/passwd entry inside the container by
# default; making these two files world-writable lets the entrypoint
# self-register the running UID/GID before handing off (same pattern as
# Red Hat's and Jupyter's "arbitrary UID"-compatible images). Safe here
# because this image runs single-tenant, locally, as whoever built or
# pulled it -- not on shared multi-tenant infrastructure.
RUN chmod 666 /etc/passwd /etc/group
COPY containers/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
