"""Snakemake profile / reproducible-workflow initialisation (`01_project_setup_and_governance/05_initialise_reproducible_workflow.py`).

Idempotent, matching the adopted scaffold-script pattern: existing files are
left untouched unless force=True. The Snakefile and config/config.yaml
already exist (created by the repository scaffold) -- this validates them
rather than re-creating them, and adds the local execution profile that
was not part of the original scaffold.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError

LOCAL_PROFILE = """cores: 8
use-conda: true
rerun-incomplete: true
printshellcmds: true
keep-going: false
latency-wait: 10
"""


def initialise_workflow(project_root: Path, force: bool = False) -> dict:
    problems = []
    for required in ["Snakefile", "config/config.yaml", "config/global_seed.yaml"]:
        if not (project_root / required).is_file():
            problems.append(required)
    if problems:
        raise PipelineError(
            f"Required workflow file(s) missing: {problems}. These are created by the repository "
            "scaffold (tools/scaffold_repository.py) or hand-authored config -- run that first."
        )

    seed_config = yaml.safe_load((project_root / "config" / "global_seed.yaml").read_text())
    if "default_seed" not in seed_config:
        raise PipelineError("config/global_seed.yaml is missing 'default_seed'.")

    profile_path = project_root / "profiles" / "local" / "config.yaml"
    wrote_profile = False
    if not profile_path.is_file() or force:
        profile_path.parent.mkdir(parents=True, exist_ok=True)
        profile_path.write_text(LOCAL_PROFILE)
        wrote_profile = True

    for d in ["results/logs", "results/tables", "results/figures"]:
        (project_root / d).mkdir(parents=True, exist_ok=True)

    return {
        "snakefile_validated": True,
        "config_validated": True,
        "default_seed": seed_config["default_seed"],
        "profile_written": wrote_profile,
        "profile_path": str(profile_path),
    }
