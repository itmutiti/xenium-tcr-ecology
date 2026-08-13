"""End-to-end reproduction smoke test (Phase 17.08).

**Scope:** this project's full pipeline (17 phases, hundreds of
scripts) took many hours of compute to run across development -- a
full end-to-end reproduction on the full dataset is out of scope for
this milestone, which instead runs a representative subset of fast,
non-external-data-dependent scripts spanning multiple phases
end-to-end (matching the scaffold's two-stage design -- "first on a
test subset" -- the "then the full dataset" stage is not attempted
here).

**Finding this process surfaced:** running this smoke test inside the
project's Docker container image (`containers/Dockerfile`) exposed a
reproducibility gap -- `environment/conda/environment.lock`
(`01_project_setup_and_governance/03_lock_software_environments.sh`'s output, last regenerated 2026-07-10, before `r-arrow`
and `r-testthat` were installed into the working environment during
later development) was missing both packages, meaning the frozen
container as originally built could not run any of this project's R
scripts or R tests. Fixed by re-running `01_project_setup_and_governance/03_lock_software_environments.sh`'s existing lock script
against the current working environment.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Representative, fast, non-external-data-dependent scripts spanning
# multiple phases -- chosen because they run in seconds, not hours, and
# do not require an external network download.
SMOKE_TEST_SCRIPTS: list[dict] = [
    {
        "phase": "14_spatial_interactions_and_barriers",
        "script": "scripts/14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py",
        "runner": "python",
    },
    {
        "phase": "15_hpv_stratified_analysis",
        "script": "scripts/15_hpv_stratified_analysis/00_validate_hpv_metadata_and_probe_signal.py",
        "runner": "python",
    },
    {
        "phase": "16_external_validation_and_generalisation",
        "script": "scripts/16_external_validation_and_generalisation/00_define_validation_claims.py",
        "runner": "python",
    },
    {
        "phase": "17_statistical_closure_and_release",
        "script": "scripts/17_statistical_closure_and_release/00_freeze_primary_results.py",
        "runner": "python",
    },
    {
        "phase": "17_statistical_closure_and_release",
        "script": "scripts/17_statistical_closure_and_release/04_generate_results_tables.py",
        "runner": "python",
    },
]


def classify_run_result(returncode: int) -> str:
    """Pure, testable: PASS for a zero exit code, FAIL otherwise
    (including a negative code from a terminating signal)."""
    return "PASS" if returncode == 0 else "FAIL"


def build_reproduction_test(project_root: Path) -> dict:
    output_dir = project_root / "reports" / "reproduction_test"
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for entry in SMOKE_TEST_SCRIPTS:
        script_path = project_root / entry["script"]
        if not script_path.is_file():
            raise PipelineError(f"'{script_path}' not found.")
        command = (
            [sys.executable, str(script_path)]
            if entry["runner"] == "python"
            else ["Rscript", str(script_path)]
        )
        result = subprocess.run(
            command, cwd=project_root, capture_output=True, text=True, timeout=300
        )
        rows.append(
            {
                "phase": entry["phase"],
                "script": entry["script"],
                "returncode": result.returncode,
                "status": classify_run_result(result.returncode),
                "stderr_tail": result.stderr[-500:] if result.returncode != 0 else "",
            }
        )

    result_df = pd.DataFrame(rows)
    output_path = output_dir / "smoke_test_results.tsv"
    result_df.to_csv(output_path, sep="\t", index=False)

    n_pass = int((result_df["status"] == "PASS").sum())

    return {
        "n_scripts_tested": len(result_df),
        "n_pass": n_pass,
        "n_fail": len(result_df) - n_pass,
        "output_path": str(output_path),
        "scope_note": (
            "Smoke test only: confirms these scripts import, execute and write their "
            "declared output on whatever derived data already exists in this checkout. "
            "Does NOT run the full DAG, does NOT use external datasets, and does NOT "
            "check any scientific result value -- see docs/execution_manual/"
            "EXECUTION_MANUAL.md's 'Smoke test vs. full reproduction' section."
        ),
    }
