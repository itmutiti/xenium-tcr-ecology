"""Project-root resolution without hardcoding a machine-specific path.

Precedence is --project-root CLI argument, then an environment variable,
then walking up from the caller's location to a marker file. This keeps
every script portable across machines and collaborators.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from xenium_tcr_ecology.infra.exceptions import PipelineError

PROJECT_ROOT_MARKER = Path("manifests/project_paths.yaml")
ENV_ROOT_VAR = "XENIUM_TCR_ECOLOGY_ROOT"


def find_project_root(cli_arg: Optional[str] = None, start: Optional[Path] = None) -> Path:
    """Resolve the repository root.

    Precedence:
      1. `cli_arg` (typically argparse's --project-root)
      2. the XENIUM_TCR_ECOLOGY_ROOT environment variable
      3. walking up from `start` (default: this file's location) looking for
         the manifests/project_paths.yaml marker file.

    Raises PipelineError, not a bare exception, so callers can catch it at
    the CLI boundary and print a clean [ERROR] message.
    """
    candidates: list[tuple[str, Path]] = []
    if cli_arg:
        candidates.append(("--project-root", Path(cli_arg)))
    env_root = os.environ.get(ENV_ROOT_VAR)
    if env_root:
        candidates.append((ENV_ROOT_VAR, Path(env_root)))

    for source, candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if (candidate / PROJECT_ROOT_MARKER).is_file():
            return candidate
        raise PipelineError(
            f"Project root supplied via {source} ('{candidate}') does not contain "
            f"'{PROJECT_ROOT_MARKER}'. Pass the path to the repository root."
        )

    here = (start or Path(__file__)).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / PROJECT_ROOT_MARKER).is_file():
            return parent

    raise PipelineError(
        "Could not locate the project root. None of the following identified it: "
        f"--project-root, ${ENV_ROOT_VAR}, or a '{PROJECT_ROOT_MARKER}' marker file "
        "found by walking up from this module. Pass --project-root explicitly or set "
        f"${ENV_ROOT_VAR}."
    )
