#!/usr/bin/env python3
"""
Verifies that every `workflow/rules/*.smk` rule's docstring-declared
"Primary output" actually exists on disk whenever that rule's sentinel has
been touched.

Why this exists: every Snakemake rule in this repository declares only a
sentinel touch-file (`results/logs/<stage>/.sentinels/NN_*.done`) as its
`output:` -- the scientific artefact itself (a parquet file, a PDF
report, a config file) is named only in the rule's docstring, not tracked
by Snakemake's own dependency graph. That is a deliberate sentinel-based
orchestration design, but it means Snakemake has no way to notice if a
script's output silently failed to get created -- only that its
sentinel got touched. This script closes that specific gap. It is not a
Snakemake rule itself; it is a post-hoc auditor, meant to be wired into
`scripts/17_statistical_closure_and_release/08_run_end_to_end_
reproduction_test.sh`'s existing "prove things actually ran" pattern.
It checks existence, not freshness -- a stale artefact from a prior run
with a still-present sentinel is not distinguished from a current one.

Handles the variety found in these docstrings (surveyed across all 139
stage rules before writing this):
  - a single literal path ("data/derived/x.parquet")
  - multiple paths separated by "; " ("data/raw/SHA256SUMS; reports/
    integrity.tsv") -- every one is checked
  - templated/glob paths with a placeholder segment ("data/staged/
    <sample_id>/") -- checked via glob, at least one match required
  - prose descriptions that are not literal paths at all ("updated h5ad/zarr
    objects; join audit") -- reported as not-checkable, never silently
    treated as passing
  - the two rules explicitly marked "NOT YET IMPLEMENTED" in their own
    docstring (the deliberately-deferred `publication_release`-only
    stubs) -- reported as deferred, never treated as failures

Usage:
    python3 tools/verify_declared_outputs.py [--project-root PATH]

Exit code 0 if every checkable, sentinel-completed rule's declared output
exists; 1 if any is missing -- a gap between "the DAG says this
finished" and "the promised artefact is actually there".
"""

from __future__ import annotations

import argparse
import glob as globmod
import re
import sys
from pathlib import Path

TOOLS_DIR = Path(__file__).resolve().parent
DEFAULT_ROOT = TOOLS_DIR.parent

RULE_RE = re.compile(
    r'rule (?P<name>\w+):\s*"""(?P<doc>.*?)"""' r'.*?output:\s*touch\("(?P<sentinel>[^"]+)"\)',
    re.DOTALL,
)
PRIMARY_OUTPUT_RE = re.compile(
    r"Primary outputs?:?\s*(?P<paths>.*?)\.(?:\s+[A-Z]|\s*$)",
    re.DOTALL,
)


def extract_rules(smk_path: Path) -> list[dict]:
    text = smk_path.read_text()
    rules = []
    for m in RULE_RE.finditer(text):
        doc = m.group("doc")
        pm = PRIMARY_OUTPUT_RE.search(doc)
        rules.append(
            {
                "file": smk_path.name,
                "name": m.group("name"),
                "sentinel": m.group("sentinel"),
                "declared_outputs": pm.group("paths").strip() if pm else None,
                "deferred": "NOT YET IMPLEMENTED" in doc,
            }
        )
    return rules


def check_output(root: Path, candidate: str) -> tuple[str, bool]:
    """Returns (kind, exists) where kind is 'literal', 'glob', or 'prose'."""
    if " " in candidate:
        return "prose", False
    if "<" in candidate:
        pattern = re.sub(r"<[^>]+>", "*", candidate).rstrip("/")
        matches = globmod.glob(str(root / pattern))
        return "glob", bool(matches)
    return "literal", (root / candidate.rstrip("/")).exists()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=str(DEFAULT_ROOT))
    args = parser.parse_args()
    root = Path(args.project_root).resolve()

    smk_files = sorted(
        p for p in (root / "workflow" / "rules").glob("*.smk") if p.name != "_checkpoints.smk"
    )

    n_pass = n_fail = n_not_run = n_not_checkable = n_deferred = 0
    failures: list[str] = []

    for smk_path in smk_files:
        for rule in extract_rules(smk_path):
            if rule["deferred"]:
                n_deferred += 1
                continue

            sentinel_path = root / rule["sentinel"]
            if not sentinel_path.exists():
                n_not_run += 1
                continue

            if not rule["declared_outputs"]:
                n_not_checkable += 1
                continue

            candidates = [c.strip() for c in rule["declared_outputs"].split(";") if c.strip()]
            results = [check_output(root, c) for c in candidates]
            checkable = [
                (c, exists) for (kind, exists), c in zip(results, candidates) if kind != "prose"
            ]

            if not checkable:
                n_not_checkable += 1
                continue

            missing = [c for c, exists in checkable if not exists]
            if missing:
                n_fail += 1
                failures.append(
                    f"  [FAIL] {rule['file']}::{rule['name']} -- sentinel done, "
                    f"but missing: {', '.join(missing)}"
                )
            else:
                n_pass += 1

    print(
        f"verify_declared_outputs: {n_pass} pass, {n_fail} fail, "
        f"{n_not_run} not-yet-run, {n_not_checkable} not-checkable (prose), "
        f"{n_deferred} deferred (NOT YET IMPLEMENTED)"
    )
    if failures:
        print("\nFailures (sentinel touched, declared output missing):")
        for line in failures:
            print(line)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
