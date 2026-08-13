"""Structured, one-JSON-file-per-invocation provenance logging.

Uses a single schema for both the per-script run log and the sentinel-JSON
provenance record written at each Snakemake checkpoint.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


def _git_commit(project_root: Path) -> Optional[str]:
    try:
        result = subprocess.run(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


class JsonRunLogger:
    """Writes one timestamped structured JSON run log per script invocation
    under results/logs/<phase_folder>/."""

    def __init__(self, logs_dir: Path, script_name: str, project_root: Path, phase: str):
        self.logs_dir = logs_dir
        self.script_name = script_name
        self.project_root = project_root
        self.phase = phase
        self.started_at = datetime.now(timezone.utc)
        self.events: list[dict] = []
        self.errors: list[str] = []
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def log_event(self, **fields: Any) -> None:
        event = dict(fields)
        event["timestamp"] = datetime.now(timezone.utc).isoformat()
        self.events.append(event)

    def log_error(self, message: str) -> None:
        self.errors.append(message)

    def write(self, status: str, extra: Optional[dict] = None) -> Path:
        record = {
            "timestamp": self.started_at.isoformat(),
            "script": self.script_name,
            "script_path": str(Path(sys.argv[0]).resolve()) if sys.argv and sys.argv[0] else None,
            "project_root": str(self.project_root),
            "phase": self.phase,
            "environment": {
                "python_version": platform.python_version(),
                "platform": platform.platform(),
                "git_commit": _git_commit(self.project_root),
            },
            "events": self.events,
            "status": status,
            "errors": self.errors,
        }
        if extra:
            record.update(extra)

        stamp = self.started_at.strftime("%Y%m%dT%H%M%SZ")
        out_path = self.logs_dir / f"{self.script_name}_{stamp}.json"
        with out_path.open("w") as fh:
            json.dump(record, fh, indent=2, default=str)
        return out_path
