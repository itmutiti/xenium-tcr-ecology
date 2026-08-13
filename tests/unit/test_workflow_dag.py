"""Regression coverage for defects surfaced by the two independent Vast.ai
clean-room verification runs: the phase01_00/phase12_04 Snakemake
orchestration flag-passthrough gaps (since resolved by removing the
underlying approval-gate concept entirely -- this project has no
committee or supervisory sign-off step), and documentation job-count
drift.

These tests shell out to a real `snakemake` binary so they exercise the
actual DAG, not a mock of it. Resolution prefers whatever `snakemake` is
already active on PATH -- inside Docker/Apptainer this is the container's
own conda env, and the documented native `.venv` is a fixed absolute path
baked into its scripts' shebang lines at creation time, which breaks under
a bind mount at a different path (e.g. Docker's `/workspace`). Falling
back to this project's documented native orchestration `.venv` (see
README.md "Setup") covers a bare native checkout where nothing is
activated on PATH yet. Skips gracefully wherever neither is present,
rather than failing on an unrelated setup gap.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from xenium_tcr_ecology.infra.paths import find_project_root

PROJECT_ROOT = find_project_root()
EXECUTION_MANUAL = PROJECT_ROOT / "docs" / "execution_manual" / "EXECUTION_MANUAL.md"


def _resolve_snakemake_bin() -> Path | None:
    on_path = shutil.which("snakemake")
    if on_path is not None:
        return Path(on_path)
    venv_bin = PROJECT_ROOT / ".venv" / "bin" / "snakemake"
    return venv_bin if venv_bin.is_file() else None


_SNAKEMAKE_BIN = _resolve_snakemake_bin()

requires_snakemake = pytest.mark.skipif(
    _SNAKEMAKE_BIN is None,
    reason=(
        "no usable snakemake found -- neither active on PATH (the case "
        "inside Docker/Apptainer) nor at "
        f"{PROJECT_ROOT / '.venv' / 'bin' / 'snakemake'} (this project's "
        "documented native orchestration .venv; see README.md 'Setup')."
    ),
)


def SNAKEMAKE_BIN() -> Path:
    """Only called from inside @requires_snakemake-guarded tests, where
    _resolve_snakemake_bin() having returned None already skipped the test."""
    assert _SNAKEMAKE_BIN is not None
    return _SNAKEMAKE_BIN


def _dry_run_total(*targets: str, forceall: bool = True) -> int:
    """Runs a real Snakemake dry run and parses the final 'total N' job count."""
    cmd = [str(SNAKEMAKE_BIN())]
    if forceall:
        cmd.append("--forceall")
    cmd += ["-n", *targets]
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, (
        f"snakemake dry run failed (exit {result.returncode}):\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    totals = re.findall(r"^total\s+(\d+)\s*$", result.stdout, re.MULTILINE)
    assert totals, f"no 'total N' line found in dry-run output:\n{result.stdout}"
    return int(totals[-1])


@requires_snakemake
class TestQuickStartCommandNeedsNoManualIntervention:
    """The literal README/EXECUTION_MANUAL quick-start command
    (`snakemake --cores N`) must resolve and, for `project_ready`
    specifically, must not require any private file or manually-supplied
    CLI flag -- this is the exact command a public reviewer runs first."""

    def test_default_dag_dry_run_resolves_with_no_manual_flags(self):
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "-n"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert (
            result.returncode == 0
        ), f"`snakemake -n` (the literal quick-start dry run) failed:\n{result.stderr}"

    def test_charter_builder_has_no_approval_gate(self):
        """Regression test: build_project_charter must remain an
        unconditional, always-succeeding record -- no acknowledge_provisional
        parameter, no committee/signoff dependency reintroduced. This project
        has no committee or supervisory approval process; nothing should ever
        require one to execute the computational workflow."""
        import inspect

        from xenium_tcr_ecology.governance.charter import build_project_charter

        params = inspect.signature(build_project_charter).parameters
        assert "acknowledge_provisional" not in params, (
            "build_project_charter has regained an acknowledge_provisional "
            "parameter -- this reintroduces an approval gate this project "
            "does not have"
        )

    def test_project_ready_dry_run_needs_no_extra_arguments(self):
        # A clean --forceall -n dry run to project_ready must succeed outright;
        # if phase01_00 ever regresses to requiring a manual pre-run, Snakemake
        # itself would still report success here (dry runs don't execute
        # shell:), so this is a structural complement to the shell: check
        # above, not a substitute for it.
        total = _dry_run_total("project_ready")
        assert total == 7, f"expected 7 jobs for project_ready, got {total}"

    def test_project_ready_does_not_require_docker_build(self):
        # Regression test: phase01_04_build_container_images hard-fails
        # without a reachable Docker daemon, confirmed absent on both
        # Vast.ai clean-room instances used to date. Requiring it here
        # would block the entire computational DAG on packaging
        # infrastructure nothing scientific depends on.
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "--forceall", "-n", "project_ready"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0
        assert "phase01_04_build_container_images" not in result.stdout, (
            "project_ready's dry-run plan still schedules phase01_04 -- this "
            "rule cannot complete without a Docker daemon and must not block "
            "the computational DAG"
        )

    def test_phase01_04_remains_independently_targetable(self):
        # The rule itself must still exist and resolve as a standalone
        # target -- it is a real, available action on hosts with genuine
        # Docker support, just not a dependency of anything else.
        total = _dry_run_total("phase01_04_build_container_images")
        assert total > 0


@requires_snakemake
class TestPhase12HasNoSeparateApprovalStep:
    """Regression test: the taxonomy-freeze decision is fully computed by
    12.03's predeclared, automated rule. There must be no separate rule
    that attributes it to a named human reviewer -- this project has no
    committee or supervisory approval process."""

    def test_external_checkpoint_resolves_with_no_manual_flags(self):
        total = _dry_run_total("external_checkpoint")
        assert total > 0

    def test_no_freeze_or_revise_decision_recording_rule_exists(self):
        result = subprocess.run(
            [
                str(SNAKEMAKE_BIN()),
                "--forceall",
                "-n",
                "phase12_04_record_freeze_or_revise_decision",
            ],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode != 0, (
            "phase12_04_record_freeze_or_revise_decision still resolves as a "
            "Snakemake target -- this project has no committee/supervisory "
            "sign-off step; the rule should not exist"
        )


@requires_snakemake
class TestCompanionReferenceAcquisitionGatesItsRealConsumers:
    """Regression test: GSE287301 (companion scRNA-seq + VDJ) previously
    had zero acquisition/verification coverage in the DAG at all --
    phase06_03 and phase08_09 read it directly with no upstream
    dependency. phase06_08 must gate both real consumers."""

    def test_phase06_08_resolves_as_a_target(self):
        total = _dry_run_total("phase06_08_acquire_companion_scrna_and_vdj_reference")
        assert total > 0

    def test_annotation_release_requires_phase06_08(self):
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "--forceall", "-n", "annotation_release"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0
        assert "phase06_08_acquire_companion_scrna_and_vdj_reference" in result.stdout

    def test_tumour_and_tcr_release_requires_phase06_08(self):
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "--forceall", "-n", "tumour_and_tcr_release"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0
        assert "phase06_08_acquire_companion_scrna_and_vdj_reference" in result.stdout
        assert "phase08_09_validate_probe_clones_against_paired_vdj_ground_truth" in result.stdout


@requires_snakemake
class TestEveryDatasetAcquisitionGatesAllItsRealConsumers:
    """Regression test for a defect *class*, not a single instance:
    a rule that reads an externally-acquired file directly, without the
    Snakemake `input:` dependency on the rule that actually acquires it,
    races under real parallel scheduling (the file may not exist, or may
    still be mid-download, when the consumer starts). This has now
    happened three times independently -- phase02's raw-archive race
    (the original Vast.ai clean-room finding), phase05_03's missing
    dependency on phase05_02 (found during the second clean-room run),
    and phase12_01's missing dependency on phase12_00 (found during the
    third, migrated-instance clean-room run: `phase12_01_test_
    transcriptional_program_transfer` reads `data/external/GSE103322/
    GSE103322_HNSCC_all_data.txt.gz` directly via `src/xenium_tcr_
    ecology/external_checkpoint/program_transfer.py`'s
    `build_program_transfer_test`, but GSE103322's acquisition is
    centralised in phase12_00 -- the only caller of `bulk_reference.
    ensure_gse103322_acquired` -- and phase12_01 never declared that
    dependency).

    Rather than adding one more one-off regression test, every
    `ensure_*_acquired`-style acquisition function in the codebase was
    audited for this same gap at the time this test was written (see
    `src/xenium_tcr_ecology/validation/spatial_dataset_acquisition.py`,
    `companion_reference_acquisition.py`, `scrna_reference_acquisition.py`,
    and `external_checkpoint/bulk_reference.py`). This table is every
    (acquisition-owner rule -> real consumer rule) edge found across all
    of them; GSE287301's two consumers already have their own dedicated
    coverage in `TestCompanionReferenceAcquisitionGatesItsRealConsumers`
    above and are intentionally not duplicated here."""

    ACQUISITION_OWNER_AND_CONSUMERS = [
        (
            "phase16_01_acquire_independent_spatial_dataset",
            "phase16_05_validate_framework_on_independent_dataset",
        ),
        (
            "phase16_08_acquire_second_independent_spatial_dataset",
            "phase16_09_validate_framework_on_second_cancer_type",
        ),
        ("phase16_02_acquire_hnscc_scrna_references", "phase16_03_validate_cell_state_signatures"),
        (
            "phase16_02_acquire_hnscc_scrna_references",
            "phase16_04_validate_ecosystem_signatures_in_bulk",
        ),
        (
            "phase12_00_project_provisional_signatures_to_bulk_reference",
            "phase12_01_test_transcriptional_program_transfer",
        ),
        (
            "phase12_00_project_provisional_signatures_to_bulk_reference",
            "phase12_02_quantify_directional_consistency",
        ),
        (
            "phase12_00_project_provisional_signatures_to_bulk_reference",
            "phase12_05_rescore_cycling_state_with_primary_method",
        ),
    ]

    @pytest.mark.parametrize("owner,consumer", ACQUISITION_OWNER_AND_CONSUMERS)
    def test_consumer_declares_owner_as_a_dependency(self, owner, consumer):
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "--forceall", "-n", consumer],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, f"dry run targeting {consumer} failed:\n{result.stderr}"
        assert owner in result.stdout, (
            f"{consumer}'s dry-run plan does not schedule {owner} -- the "
            "declared Snakemake `input:` dependency on the rule that "
            "acquires this consumer's external dataset appears to be "
            "missing or has regressed, which will race under real "
            "parallel scheduling exactly as phase12_01 did"
        )


@requires_snakemake
class TestApptainerIsOptionalAndDoesNotAlterTheNativeDag:
    """Regression test: adding Apptainer as an optional execution route
    (containers/Apptainer.def, profiles/apptainer/, tools/run_with_
    apptainer.sh) must not change the native scientific rule graph, add
    a container-build rule to the default DAG, or require Apptainer to
    be installed for native execution."""

    def test_native_forceall_job_count_unchanged(self):
        # The single `container: "xenium-tcr-ecology.sif"` directive
        # added to the Snakefile must be completely inert for native
        # execution -- same job count as before Apptainer was added.
        assert _dry_run_total() == 152

    def test_no_container_build_rule_in_default_dag(self):
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "--forceall", "-n"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0
        lowered = result.stdout.lower()
        assert "apptainer" not in lowered and "sif" not in lowered, (
            "the default DAG must not reference building or fetching a "
            "container image -- Apptainer is an optional execution route, "
            "not a prerequisite for the computational workflow"
        )

    def test_native_dry_run_succeeds_regardless_of_apptainer_availability(self):
        # This test environment may or may not have `apptainer` on PATH;
        # either way, native dry-run resolution must succeed unmodified.
        result = subprocess.run(
            [str(SNAKEMAKE_BIN()), "-n"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0

    def test_apptainer_profile_is_not_the_default(self):
        # Confirms the profile is a separate, explicitly-opted-into
        # config, never picked up by a bare `snakemake` invocation.
        profile_path = PROJECT_ROOT / "profiles" / "apptainer" / "config.v8+.yaml"
        assert profile_path.is_file()
        default_profile_dir = PROJECT_ROOT / ".config" / "snakemake"
        assert (
            not default_profile_dir.exists() or not (default_profile_dir / "config.yaml").exists()
        ), "no config should make the apptainer profile Snakemake's implicit default"

    def test_launcher_script_is_a_separate_opt_in_entry_point(self):
        launcher = PROJECT_ROOT / "tools" / "run_with_apptainer.sh"
        assert launcher.is_file()
        # Never *invoked* by any rule's shell: command -- a documentation
        # cross-reference (e.g. in the Snakefile's own module docstring)
        # is fine; an actual `shell:` line calling it would not be.
        shell_line_pattern = re.compile(r'shell:\s*".*run_with_apptainer')
        for smk in (PROJECT_ROOT / "workflow" / "rules").glob("*.smk"):
            assert not shell_line_pattern.search(smk.read_text())


@requires_snakemake
class TestDocumentedJobCountMatchesLiveDag:
    """Regression test for the off-by-one job-count drift the forensic audit
    found: EXECUTION_MANUAL.md's cited dry-run job count must match a real,
    from-clean-state `snakemake --forceall -n` total, not merely whatever
    number was true when the sentence was last hand-edited."""

    def test_execution_manual_job_count_matches_live_dag(self):
        assert EXECUTION_MANUAL.is_file()
        text = EXECUTION_MANUAL.read_text()
        match = re.search(r"resolves to the same\s*(\d+)\s*jobs", text)
        assert match, (
            "could not find the '...resolves to the same N jobs...' sentence in "
            "EXECUTION_MANUAL.md -- has its wording changed? Update this test's "
            "pattern to match, not just the number."
        )
        documented_count = int(match.group(1))
        live_count = _dry_run_total()
        assert documented_count == live_count, (
            f"EXECUTION_MANUAL.md documents {documented_count} jobs for a full "
            f"`--forceall -n` dry run, but the live DAG currently resolves to "
            f"{live_count}. Update the manual's figure (and regenerate its PDF)."
        )
