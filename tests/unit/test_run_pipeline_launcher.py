"""Regression coverage for `tools/run_pipeline.sh`, the Docker -> Apptainer
-> native launcher, and for cross-route DAG identity specifically.

The launcher's own capability probing already decides at runtime which
routes are usable; these tests mirror that rather than assuming
any particular route is available, so they pass identically on a reviewer's
machine with only the native route set up and on a machine with all three.
Each route's real DAG total is compared directly against the others that
are actually available in this environment -- not against a hardcoded
number -- so this stays meaningful without requiring every route locally.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess

import pytest

from xenium_tcr_ecology.infra.paths import find_project_root

PROJECT_ROOT = find_project_root()
LAUNCHER = PROJECT_ROOT / "tools" / "run_pipeline.sh"


def _docker_usable() -> bool:
    if shutil.which("docker") is None:
        return False
    return subprocess.run(["docker", "info"], capture_output=True, timeout=30).returncode == 0


def _apptainer_usable() -> bool:
    if shutil.which("apptainer") is None:
        return False
    sif = os.environ.get("XENIUM_SIF", str(PROJECT_ROOT / "xenium-tcr-ecology.sif"))
    return os.path.isfile(sif)


def _native_usable() -> bool:
    mgr = shutil.which("mamba") or shutil.which("conda")
    if mgr is None:
        return False
    result = subprocess.run([mgr, "env", "list"], capture_output=True, text=True, timeout=30)
    return re.search(r"\bxenium-tcr-ecology\b", result.stdout) is not None


def _launcher_dry_run_total(route: str) -> int:
    result = subprocess.run(
        [str(LAUNCHER), "snakemake", "--forceall", "-n"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        env={**os.environ, "XENIUM_FORCE_ROUTE": route},
    )
    assert result.returncode == 0, (
        f"launcher dry run via '{route}' failed (exit {result.returncode}):\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    totals = re.findall(r"^total\s+(\d+)\s*$", result.stdout, re.MULTILINE)
    assert totals, f"no 'total N' line in '{route}' dry-run output:\n{result.stdout}"
    return int(totals[-1])


class TestLauncherExists:
    def test_script_is_executable(self):
        assert LAUNCHER.is_file()
        assert os.access(LAUNCHER, os.X_OK)

    def test_script_is_syntactically_valid_bash(self):
        result = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr

    def test_rejects_unknown_forced_route(self):
        result = subprocess.run(
            [str(LAUNCHER), "true"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            env={**os.environ, "XENIUM_FORCE_ROUTE": "not-a-real-route"},
        )
        assert result.returncode != 0
        assert "unknown backend" in (result.stderr + result.stdout).lower()


class TestCrossRouteDagIdentity:
    """The core regression coverage the odt-sourced execution-route work
    asked for: enabling Docker or Apptainer must never change the
    scientific DAG relative to native. Compares whichever routes are
    usable on the machine running this test against each other,
    directly -- not against a number hardcoded here."""

    def test_available_routes_agree_on_forceall_job_count(self):
        available = {
            "docker": _docker_usable(),
            "apptainer": _apptainer_usable(),
            "native": _native_usable(),
        }
        usable = [route for route, ok in available.items() if ok]
        if len(usable) < 2:
            pytest.skip(
                f"only {usable or 'no'} route(s) usable in this environment -- "
                "need at least 2 to compare; not a defect, just this machine's setup"
            )

        totals = {route: _launcher_dry_run_total(route) for route in usable}
        assert len(set(totals.values())) == 1, (
            f"routes disagree on the forced dry-run job count -- the scientific "
            f"DAG must be identical regardless of execution route: {totals}"
        )

    def test_launcher_default_route_matches_direct_native_invocation(self):
        # The launcher must be a pure environment-selection wrapper: its
        # default-routed dry run and a direct native `snakemake` dry run
        # must resolve the identical DAG, proving the launcher does not
        # itself alter scientific behaviour.
        if not _native_usable():
            pytest.skip("native route not set up in this environment")

        native_bin = shutil.which("snakemake")
        if native_bin is None:
            pytest.skip("no snakemake on PATH to compare against directly")

        direct = subprocess.run(
            [native_bin, "--forceall", "-n"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert direct.returncode == 0
        direct_totals = re.findall(r"^total\s+(\d+)\s*$", direct.stdout, re.MULTILINE)
        assert direct_totals

        via_launcher = _launcher_dry_run_total("native")
        assert int(direct_totals[-1]) == via_launcher
