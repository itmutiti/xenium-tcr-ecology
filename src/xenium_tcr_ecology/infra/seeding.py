"""Central, documented random-seed configuration.

Every stochastic procedure in this project (permutation tests, bootstrap
confidence intervals, subsampling, k-means initialisation) was already
seeded -- verified by re-running each
identified stochastic script twice and comparing checksums: every one
produced byte-identical output. What was missing was centralisation: the
same seed values (`42` and, in one historical sub-lineage, `0`) were
hardcoded as a duplicated local constant (`RNG_SEED = 42`, `rng_seed:
int = 0`, ...) independently in ~30 files, rather than read from one
source of truth. This module fixes that, without changing any of those
already-frozen, checksummed numeric values.

Two distinct seeds exist in this codebase's history, both preserved
exactly as originally used rather than silently unified:

- `default_seed` (42) -- the general-purpose seed used by every
  permutation test, bootstrap CI, and subsampling procedure outside the
  annotation stack (graphs, clone ecology, interactions, checkpoint,
  HPV, validation modules).
- `annotation_seed` (0) -- used consistently throughout the annotation/
  preprocessing/tumour-scoring stack (stages 05-07: clustering, T-cell
  substate scoring, program scores, malignancy scoring, blinded review),
  predating `default_seed`'s introduction elsewhere. Kept as its own
  named constant rather than retroactively changed to 42: doing so would
  silently alter already-computed, checksummed, deeply downstream-
  dependent results (most of stages 08-17 depend on stage 06's
  annotation output) for no reproducibility benefit -- 0 is exactly as
  reproducible as 42, as long as it is used consistently, which this
  module now guarantees structurally rather than by convention.

Both are read from `config/config.yaml`, cached per-process (config is
read once, not once per call).
"""

from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path

import yaml

from xenium_tcr_ecology.infra.paths import find_project_root


@lru_cache(maxsize=1)
def _load_seed_config() -> dict:
    project_root = find_project_root(start=Path(__file__))
    config_path = project_root / "config" / "config.yaml"
    return yaml.safe_load(config_path.read_text())


def get_default_seed() -> int:
    """The general-purpose seed (42) -- use for any new stochastic
    procedure outside the annotation stack."""
    return int(_load_seed_config()["default_seed"])


def get_annotation_seed() -> int:
    """The annotation-stack seed (0) -- use only when matching the
    existing, already-frozen annotation/preprocessing/tumour-scoring
    stack's historical seed; see module docstring."""
    return int(_load_seed_config()["annotation_seed"])


def deterministic_seed(*parts: object) -> int:
    """Process-stable, deterministic seed derived from arbitrary values
    (patient id, effect size, ...) -- unlike Python's built-in `hash()`
    on a tuple containing any `str` element (randomised per process by
    default), `hashlib.md5` on a stable string encoding of `parts`
    gives the same seed every time, in every process. Canonical
    implementation, moved here unchanged (same hash algorithm, same
    seed values) from `xenium_tcr_ecology.graphs.null_models`, which
    now re-exports this exact function rather than duplicating it."""
    key = "|".join(str(p) for p in parts).encode()
    return int(hashlib.md5(key).hexdigest(), 16) % (2**32)
