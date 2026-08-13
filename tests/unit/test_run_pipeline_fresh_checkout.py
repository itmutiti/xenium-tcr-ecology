"""Regression coverage for tools/run_pipeline.sh's backend-selection
architecture on a fresh repository checkout: no
environment/conda/environment.lock, no local Docker image, no .sif.

These tests run the real script as a subprocess (not a reimplementation
of its logic), against fake `docker`/`apptainer`/`mamba` stub binaries on
a controlled PATH, so they exercise the actual selection branches without
needing a real multi-minute image build per test. The stubs are
deliberately minimal: just enough to prove which backend the script
chose and why, not to simulate a real build's output.
"""

from __future__ import annotations

import os
import pathlib
import shutil
import stat
import subprocess

import pytest

from xenium_tcr_ecology.infra.paths import find_project_root

PROJECT_ROOT = find_project_root()

FAKE_DOCKER = """#!/bin/bash
case "$1" in
  info) exit 0 ;;
  run)
    if [[ "$*" == *"hello-world"* ]]; then exit 0; fi
    echo "FAKE_DOCKER_RAN"
    exit 0
    ;;
  image) exit 1 ;;  # always "not found" -- forces the build path every time
  build)
    echo "FAKE_DOCKER_BUILD"
    exit 0
    ;;
  *) exit 1 ;;
esac
"""

FAKE_APPTAINER = """#!/bin/bash
case "$1" in
  build)
    echo "FAKE_APPTAINER_BUILD"
    # create the .sif path (last arg before the .def) so a subsequent
    # `sha256sum` on it in record_route succeeds
    touch "${@: -2:1}" 2>/dev/null || true
    exit 0
    ;;
  exec)
    echo "FAKE_APPTAINER_RAN"
    exit 0
    ;;
  *) exit 1 ;;
esac
"""

FAKE_MAMBA = """#!/bin/bash
case "$1" in
  env)
    if [[ "$2" == "list" ]]; then
      echo "  xenium-tcr-ecology   /fake/envs/xenium-tcr-ecology"
      exit 0
    fi
    exit 0
    ;;
  run)
    echo "FAKE_MAMBA_RAN"
    exit 0
    ;;
  *) exit 1 ;;
esac
"""


_REAL_BACKEND_BINARY_NAMES = {"docker", "apptainer", "mamba", "conda"}


def _sanitized_path_dir(tmp_path):
    """Symlinks every real binary on the current PATH into one directory,
    except docker/apptainer/mamba/conda -- so tests can put this ahead of
    (and *instead of*, not alongside) the real PATH and get a host where
    those specific tools are absent, not merely shadowed, while
    bash/coreutils/python/etc. all still work normally."""
    sanitized = tmp_path / "sanitized_real_bin"
    sanitized.mkdir(exist_ok=True)
    seen = set()
    for directory in os.environ["PATH"].split(os.pathsep):
        d = pathlib.Path(directory)
        if not d.is_dir():
            continue
        for entry in d.iterdir():
            if entry.name in _REAL_BACKEND_BINARY_NAMES or entry.name in seen:
                continue
            if not os.access(entry, os.X_OK):
                continue
            seen.add(entry.name)
            try:
                (sanitized / entry.name).symlink_to(entry)
            except OSError:
                pass
    return sanitized


def _route_log(project_root):
    # docker build's/apptainer build's own output is redirected to this
    # log file (>>"$ROUTE_LOG" 2>&1), not to the launcher's own stdout --
    # only the info()/error() lines are tee'd to both.
    log = project_root / "results" / "logs" / "run_pipeline" / "route_selection.log"
    return log.read_text() if log.is_file() else ""


def _write_stub(bin_dir, name, content):
    path = bin_dir / name
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _run_launcher(fake_bin_dir, sanitized_dir, project_root, *args):
    # Invoke *this copy's own* tools/run_pipeline.sh, not the real
    # repository's -- the script resolves its project root from its own
    # BASH_SOURCE[0], not from cwd, so running the real script with cwd
    # set to a fake checkout would silently still operate on the real
    # repository and defeat the whole point of this fixture.
    launcher = project_root / "tools" / "run_pipeline.sh"
    # PATH is built *only* from the fake stubs and the sanitized
    # (docker/apptainer/mamba/conda-excluded) real-binary dir -- NOT the
    # raw real PATH -- otherwise the real docker/apptainer/mamba actually
    # installed on this machine would still be reachable regardless of
    # what is or isn't stubbed, defeating these tests' whole purpose.
    env = {**os.environ, "PATH": os.pathsep.join([str(fake_bin_dir), str(sanitized_dir)])}
    return subprocess.run(
        [str(launcher), *args],
        cwd=project_root,
        capture_output=True,
        text=True,
        timeout=30,
        env=env,
    )


@pytest.fixture
def fresh_checkout(tmp_path):
    """A real copy of the repository's tracked files (via `git ls-files`
    -- respects .gitignore automatically) into a scratch directory, with
    environment.lock, any .sif, and any snakemake state removed, so it
    has none of the generated artefacts that used to gate
    backend selection."""
    dest = tmp_path / "fresh_checkout"
    dest.mkdir()
    tracked = subprocess.run(
        ["git", "ls-files"], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True
    ).stdout.splitlines()
    for rel in tracked:
        src = PROJECT_ROOT / rel
        if not src.is_file():
            continue
        dst = dest / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    lock = dest / "environment" / "conda" / "environment.lock"
    if lock.exists():
        lock.unlink()
    sif = dest / "xenium-tcr-ecology.sif"
    if sif.exists():
        sif.unlink()
    return dest


@pytest.fixture
def fake_bin(tmp_path):
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    return bin_dir


@pytest.fixture
def sanitized_bin(tmp_path):
    """Every real binary on PATH except docker/apptainer/mamba/conda --
    so a test can control exactly which of those four are "available" via
    fake_bin stubs, without the real ones (wherever they happen to live
    on this machine) leaking through and defeating the test."""
    return _sanitized_path_dir(tmp_path)


class TestFreshCheckoutNeverGatesOnGeneratedFiles:
    def test_docker_selected_despite_missing_environment_lock(
        self, fresh_checkout, fake_bin, sanitized_bin
    ):
        _write_stub(fake_bin, "docker", FAKE_DOCKER)
        assert not (fresh_checkout / "environment" / "conda" / "environment.lock").exists()

        result = _run_launcher(fake_bin, sanitized_bin, fresh_checkout, "true")

        assert "Backend selected: docker" in result.stdout, result.stdout + result.stderr
        assert "no environment/conda/environment.lock" not in result.stdout.lower()
        assert "FAKE_DOCKER_BUILD" in _route_log(fresh_checkout)

    def test_apptainer_selected_automatically_when_docker_unavailable(
        self, fresh_checkout, fake_bin, sanitized_bin
    ):
        # No docker stub at all, and docker excluded from sanitized_bin
        # too -- absent, not just failing.
        _write_stub(fake_bin, "apptainer", FAKE_APPTAINER)

        result = _run_launcher(fake_bin, sanitized_bin, fresh_checkout, "true")

        assert "Docker: rejected" in result.stdout, result.stdout + result.stderr
        assert "Backend selected: apptainer" in result.stdout
        assert "FAKE_APPTAINER_BUILD" in _route_log(fresh_checkout)

    def test_native_selected_only_when_both_containers_unavailable(
        self, fresh_checkout, fake_bin, sanitized_bin
    ):
        # Neither docker nor apptainer stubbed or reachable at all.
        _write_stub(fake_bin, "mamba", FAKE_MAMBA)

        result = _run_launcher(fake_bin, sanitized_bin, fresh_checkout, "true")

        assert "Docker: rejected" in result.stdout, result.stdout + result.stderr
        assert "Apptainer: rejected" in result.stdout
        assert "Backend selected: native" in result.stdout

    def test_every_rejection_and_the_final_selection_are_reported(
        self, fresh_checkout, fake_bin, sanitized_bin
    ):
        _write_stub(fake_bin, "mamba", FAKE_MAMBA)

        result = _run_launcher(fake_bin, sanitized_bin, fresh_checkout, "true")

        assert "docker client not found on PATH" in result.stdout
        assert "apptainer not found on PATH" in result.stdout
        assert "Backend selected: native" in result.stdout


class TestExplicitBackendNeverSilentlyFallsBack:
    def test_explicit_docker_request_fails_loudly_when_unavailable(
        self, fresh_checkout, fake_bin, sanitized_bin
    ):
        _write_stub(fake_bin, "mamba", FAKE_MAMBA)  # native IS available as an alternative

        result = _run_launcher(
            fake_bin, sanitized_bin, fresh_checkout, "--backend", "docker", "true"
        )

        assert result.returncode != 0
        assert "no automatic fallback has been performed" in result.stdout.lower()
        assert "Backend selected: native" not in result.stdout
        assert "native (available)" in result.stdout

    def test_explicit_native_request_fails_loudly_when_unavailable(
        self, fresh_checkout, fake_bin, sanitized_bin
    ):
        _write_stub(fake_bin, "docker", FAKE_DOCKER)  # docker IS available as an alternative
        # Deliberately no mamba/conda stub -- native unusable,
        # not just failing (sanitized_bin excludes any real conda/mamba).

        result = _run_launcher(
            fake_bin, sanitized_bin, fresh_checkout, "--backend", "native", "true"
        )

        assert result.returncode != 0
        assert "no automatic fallback has been performed" in result.stdout.lower()
        assert "docker (available)" in result.stdout

    def test_unknown_backend_rejected_immediately(self, fresh_checkout, fake_bin, sanitized_bin):
        result = _run_launcher(
            fake_bin, sanitized_bin, fresh_checkout, "--backend", "not-a-real-backend", "true"
        )
        assert result.returncode != 0
        assert "unknown backend" in (result.stdout + result.stderr).lower()
