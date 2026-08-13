"""Typed exception for conditions a script should abort on with a clear,
actionable message -- as opposed to an uncaught bug, which should surface as
a normal traceback. Every phase script catches PipelineError at its CLI
boundary and prints "[ERROR] <message>" rather than a raw stack trace."""

from __future__ import annotations


class PipelineError(RuntimeError):
    """Raised for any condition that should abort a script gracefully:
    a missing/malformed input, a failed contract validation, a policy
    violation, or a governance gate that has not been satisfied (e.g. the
    External Checkpoint Validation taxonomy-freeze gate not yet recording a 'frozen' decision)."""
