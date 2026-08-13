"""Spatial null models and calibration suite (`09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`).

Implements and calibrates the three null models named in the registered
`null_model_calibration_suite` analysis (`governance/analysis_registry.tsv`)
against `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s synthetic ground-truth data, before any of them is
trusted on real clones -- this gates all downstream spatial-association
claims (Phases 10-14), per the registry's multiplicity-family note.

Test statistic (all three null models): mean distance from each "clone"
cell to its nearest "tumour" cell -- a standard, interpretable
spatial-attraction statistic, directly analogous to `07_tumour_epithelium_characterisation/05_extract_tumour_boundaries.py`'s
signed-distance-to-boundary framework. A one-tailed permutation p-value
(fraction of null-permutation statistics <= the observed statistic) tests
for attraction (clone cells closer to tumour than expected by chance).

**Three distinct null models, a methodological judgment call
made and documented since the specification names them only
in one brief phrase:**
1. **Constrained permutation** (`free_permutation_pvalue`): clone labels
   are randomly reassigned among all eligible (non-tumour) points,
   unconstrained beyond keeping the same clone count -- the classical,
   simplest spatial permutation null.
2. **Degree-preserving** (`degree_matched_permutation_pvalue`): clone
   labels are randomly reassigned only among points with matching local
   density (their degree in `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`'s calibrated 30um graph, binned
   into quartiles) -- a degree-preserving permutation in the network-
   null-model sense: preserves the clone class's local-density profile,
   not just its count, so a spurious "attraction" driven purely by clone
   cells sitting in denser regions (not truly tumour-proximal ones) is
   correctly not rewarded.
3. **Graph-preserving** (`graph_distance_permutation_pvalue`): the same
   unconstrained label permutation as (1), but the test statistic itself
   is computed using graph-hop distance (BFS shortest path over the fixed
   graph topology) to the nearest tumour cell, not raw Euclidean distance
   -- tests whether a finding survives switching from a metric to a
   topological notion of "distance," holding the graph itself fixed.

**Unit of analysis:** one synthetic patient replicate (`09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py`'s
`unit_of_analysis: "synthetic dataset replicate"`) -- each of the 10
synthetic patients at a given effect size yields its own p-value; Type I
error and power are the empirical rejection rate across these 10
independent replicates, at nominal alpha=0.05. Explicitly caveated: 10
replicates gives a wide confidence interval on any empirical rate
estimate (a Clopper-Pearson 95% CI for a true rate of 0.05 at n=10
comfortably spans 0 to ~0.25) -- reported as such, not overclaimed as a
precise calibration.

**Reproducibility bug found and fixed by `17_statistical_closure_and_release/09_run_null_model_calibration_regression.py`'s
calibration-regression check:** `build_null_model_calibration` used to
seed each replicate's RNG via Python's built-in `hash((patient_id,
effect_size, "calibration")) % (2**32)`. Python randomises the hash of
`str` objects per process by default (`PYTHONHASHSEED`) -- confirmed,
not assumed: `hash(("SYN00", 0.0, "calibration"))` gives a
different value on every separate Python invocation. Any tuple
containing a string element is affected; a tuple of only `int`/`float`
elements is not (confirmed: `synthetic_patterns.py`'s
`hash((patient_idx, effect_size))` seeding is already deterministic
and did not need this fix). `deterministic_seed` below uses `hashlib.
md5` instead -- stable across process invocations by construction,
fixing a genuine (if statistically minor) reproducibility gap in
already-established `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py` production code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.spatial import cKDTree

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import deterministic_seed  # noqa: F401 (re-exported)
from xenium_tcr_ecology.graphs.candidate_graphs import build_radius_graph

N_PERMUTATIONS = 199
NOMINAL_ALPHA = 0.05
N_DEGREE_STRATA = 4


def compute_mean_distance_to_tumour(
    points: np.ndarray, is_tumour: np.ndarray, is_clone: np.ndarray
) -> float:
    tumour_tree = cKDTree(points[is_tumour])
    distances, _ = tumour_tree.query(points[is_clone])
    return float(distances.mean())


def free_permutation_pvalue(
    points: np.ndarray,
    is_tumour: np.ndarray,
    is_clone: np.ndarray,
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
) -> float:
    observed = compute_mean_distance_to_tumour(points, is_tumour, is_clone)
    eligible = np.where(~is_tumour)[0]
    n_clone = int(is_clone.sum())
    null_stats = np.empty(n_permutations)
    for i in range(n_permutations):
        chosen = rng.choice(eligible, size=n_clone, replace=False)
        perm_is_clone = np.zeros(len(points), dtype=bool)
        perm_is_clone[chosen] = True
        null_stats[i] = compute_mean_distance_to_tumour(points, is_tumour, perm_is_clone)
    return float((np.sum(null_stats <= observed) + 1) / (n_permutations + 1))


def compute_degree_strata(
    points: np.ndarray, radius_um: float, n_strata: int = N_DEGREE_STRATA
) -> np.ndarray:
    """Pure, testable degree computation + quartile stratification for the
    degree-preserving null model."""
    graph = build_radius_graph(points, radius_um)
    degree = np.asarray(graph.sum(axis=1)).ravel()
    return pd.qcut(degree, q=n_strata, labels=False, duplicates="drop")


def degree_matched_permutation_pvalue(
    points: np.ndarray,
    is_tumour: np.ndarray,
    is_clone: np.ndarray,
    degree_strata: np.ndarray,
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
) -> float:
    observed = compute_mean_distance_to_tumour(points, is_tumour, is_clone)
    eligible = np.where(~is_tumour)[0]
    clone_strata_counts = pd.Series(degree_strata[is_clone]).value_counts()

    null_stats = np.empty(n_permutations)
    for i in range(n_permutations):
        perm_is_clone = np.zeros(len(points), dtype=bool)
        for stratum, n_needed in clone_strata_counts.items():
            stratum_pool = eligible[degree_strata[eligible] == stratum]
            n_take = min(n_needed, len(stratum_pool))
            chosen = rng.choice(stratum_pool, size=n_take, replace=False)
            perm_is_clone[chosen] = True
        null_stats[i] = compute_mean_distance_to_tumour(points, is_tumour, perm_is_clone)
    return float((np.sum(null_stats <= observed) + 1) / (n_permutations + 1))


def graph_distance_permutation_pvalue(
    graph: sparse.csr_matrix,
    is_tumour: np.ndarray,
    is_clone: np.ndarray,
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
) -> float:
    graph_distances = shortest_path(
        graph, directed=False, unweighted=False, indices=np.where(is_tumour)[0]
    )
    min_dist_to_tumour = graph_distances.min(axis=0)
    finite = np.isfinite(min_dist_to_tumour)

    def _stat(clone_mask: np.ndarray) -> float:
        relevant = clone_mask & finite
        if not relevant.any():
            return np.inf
        return float(min_dist_to_tumour[relevant].mean())

    observed = _stat(is_clone)
    eligible = np.where(~is_tumour & finite)[0]
    n_clone = int(is_clone.sum())
    null_stats = np.empty(n_permutations)
    for i in range(n_permutations):
        chosen = rng.choice(eligible, size=min(n_clone, len(eligible)), replace=False)
        perm_is_clone = np.zeros(len(is_tumour), dtype=bool)
        perm_is_clone[chosen] = True
        null_stats[i] = _stat(perm_is_clone)
    return float((np.sum(null_stats <= observed) + 1) / (n_permutations + 1))


# `deterministic_seed` moved to xenium_tcr_ecology.infra.seeding (the
# single, central, documented seed source for the whole project) and
# re-exported here unchanged (same hashlib.md5 implementation, same
# seed values) so existing callers of
# `xenium_tcr_ecology.graphs.null_models.deterministic_seed` are unaffected.


def clopper_pearson_ci(n_successes: int, n_trials: int, alpha: float = 0.05) -> tuple[float, float]:
    from scipy.stats import beta

    if n_trials == 0:
        return (np.nan, np.nan)
    lower = (
        0.0 if n_successes == 0 else beta.ppf(alpha / 2, n_successes, n_trials - n_successes + 1)
    )
    upper = (
        1.0
        if n_successes == n_trials
        else beta.ppf(1 - alpha / 2, n_successes + 1, n_trials - n_successes)
    )
    return (float(lower), float(upper))


def build_null_model_calibration(project_root: Path) -> dict:
    synthetic_dir = project_root / "data" / "graphs" / "synthetic"
    graph_params_path = project_root / "config" / "graph_parameters.yaml"
    output_path = project_root / "reports" / "graphs" / "null_model_calibration.parquet"

    if not synthetic_dir.is_dir():
        raise PipelineError(
            f"'{synthetic_dir}' not found. Run `09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py` first."
        )
    if not graph_params_path.is_file():
        raise PipelineError(
            f"'{graph_params_path}' not found. Run `09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py` first."
        )

    radius_um = yaml.safe_load(graph_params_path.read_text())["calibrated_radius_um"]
    manifest = pd.read_csv(synthetic_dir / "_manifest.tsv", sep="\t")

    rows = []
    for _, row in manifest.iterrows():
        replicate = pd.read_parquet(row["output_path"])
        points = replicate[["x", "y"]].to_numpy()
        is_tumour = replicate["is_tumour"].to_numpy()
        is_clone = replicate["is_clone"].to_numpy()
        rng = np.random.default_rng(
            deterministic_seed(row["patient_id"], row["effect_size"], "calibration")
        )

        graph = build_radius_graph(points, radius_um)
        degree_strata = compute_degree_strata(points, radius_um)

        p_free = free_permutation_pvalue(points, is_tumour, is_clone, rng)
        p_degree = degree_matched_permutation_pvalue(
            points, is_tumour, is_clone, degree_strata, rng
        )
        p_graph = graph_distance_permutation_pvalue(graph, is_tumour, is_clone, rng)

        rows.append(
            {
                "patient_id": row["patient_id"],
                "effect_size": row["effect_size"],
                "pvalue_constrained_permutation": p_free,
                "pvalue_degree_preserving": p_degree,
                "pvalue_graph_preserving": p_graph,
            }
        )

    results = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_parquet(output_path)

    null_model_cols = [
        "pvalue_constrained_permutation",
        "pvalue_degree_preserving",
        "pvalue_graph_preserving",
    ]
    calibration_summary = {}
    for col in null_model_cols:
        by_effect = {}
        for effect_size, group in results.groupby("effect_size"):
            n_rejected = int((group[col] < NOMINAL_ALPHA).sum())
            n_total = len(group)
            ci = clopper_pearson_ci(n_rejected, n_total)
            by_effect[str(effect_size)] = {
                "rejection_rate": round(n_rejected / n_total, 4),
                "n_rejected": n_rejected,
                "n_total": n_total,
                "ci_95": [round(ci[0], 4), round(ci[1], 4)],
            }
        calibration_summary[col] = by_effect

    return {
        "n_replicates": len(results),
        "nominal_alpha": NOMINAL_ALPHA,
        "n_permutations": N_PERMUTATIONS,
        "calibration_summary": calibration_summary,
        "output_path": str(output_path),
    }
