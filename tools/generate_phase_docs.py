#!/usr/bin/env python3
"""
Generate docs/phases.md from manifests/phase_registry.yaml and the
scripts/<phase>/ directory listing, so the documentation site cannot
silently drift out of sync with the phase structure the way a hand-
maintained page would (see docs/index.md).

Only git-tracked scripts are listed: a script present on disk but
deliberately excluded from git (e.g. via `.git/info/exclude`, the same
mechanism used for other private, local-only files) must not leak into
this public-facing reference just because it still exists locally.

Usage:
    python3 tools/generate_phase_docs.py
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def _tracked_files(root: Path) -> set[str]:
    result = subprocess.run(
        ["git", "ls-files"], cwd=root, capture_output=True, text=True, check=True
    )
    return set(result.stdout.splitlines())


def main() -> int:
    registry = yaml.safe_load((ROOT / "manifests" / "phase_registry.yaml").read_text())
    tracked = _tracked_files(ROOT)

    lines = [
        "# Phase reference",
        "",
        "Generated from `manifests/phase_registry.yaml` -- do not hand-edit.",
        "",
    ]
    for folder in sorted(registry):
        entry = registry[folder]
        script_dir = ROOT / "scripts" / folder
        scripts = (
            sorted(
                p.name
                for p in script_dir.glob("*")
                if p.is_file()
                and f"scripts/{folder}/{p.name}" in tracked
            )
            if script_dir.is_dir()
            else []
        )
        lines.append(f"## `{folder}`")
        lines.append("")
        lines.append(entry["purpose"])
        lines.append("")
        lines.append(f"Deterministic: `{entry['deterministic']}`")
        lines.append("")
        if scripts:
            lines.append("Scripts:")
            lines.append("")
            for s in scripts:
                lines.append(f"- `{s}`")
            lines.append("")

    out_path = ROOT / "docs" / "phases.md"
    out_path.write_text("\n".join(lines))
    print(f"[OK] Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
