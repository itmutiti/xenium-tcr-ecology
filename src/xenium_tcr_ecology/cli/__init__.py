"""Shared CLI conventions for phase scripts.

Every phase script's module docstring doubles as its --help text via
RawDescriptionHelpFormatter, every script exposes --project-root, and
main() returns an exit code.
"""

from __future__ import annotations

import argparse


def base_parser(doc: str) -> argparse.ArgumentParser:
    """Build an ArgumentParser pre-configured with this project's house
    style: the calling script's module docstring as help text, and the
    universal --project-root flag. Phase scripts add their own arguments to
    the returned parser before calling parse_args()."""
    parser = argparse.ArgumentParser(
        description=doc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Repository root. Defaults to $XENIUM_TCR_ECOLOGY_ROOT or a marker-file walk-up.",
    )
    return parser
