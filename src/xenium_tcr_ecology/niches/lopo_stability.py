"""Leave-one-patient-out (LOPO) niche stability (`10_niche_and_ecosystem_discovery/07_leave_one_patient_out_niche_stability.py`).

Tests whether `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s archetype structure remains identifiable
when each patient is withheld entirely from the clustering fit -- the
strongest generalisation check available for the Niche and Ecosystem Discovery pipeline.
If archetypes only reflected one patient's idiosyncratic tissue, removing
that patient and refitting would produce very different centroids, and
the withheld patient's cells would be poorly predicted by the remaining
patients' structure. No data leakage risk: different patients' tissue
sections are already spatially independent (`09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s primary
graph never connects cells across sections), so a section's composition
vectors (`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`) are computed entirely from that section's
neighbours, unaffected by which other patient is held out.

Reuses `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s established K=6 (selected via the full
consensus-clustering PAC procedure) for every LOPO fold, rather than
re-running full K-selection 11 times -- computationally prohibitive, and
would risk each fold choosing a different K, breaking cross-fold
archetype correspondence. Each fold refits only the final k-means
assignment step on the n-1 remaining patients' radius_30.0um
composition vectors (same scale as `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`).

Cluster labels are arbitrary across independent k-means fits, so each
LOPO fit's centroids are matched to `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s full-dataset centroids
via the Hungarian algorithm (`scipy.optimize.linear_sum_assignment`) on
pairwise Euclidean centroid distance -- a standard technique for
cross-fit cluster correspondence, not an assumption that cluster index
i means the same archetype across fits.

Two stability metrics per fold:
1. Centroid stability: cosine similarity between each LOPO centroid and
   its matched full-dataset centroid -- do the archetypes themselves
   look the same without this patient's data?
2. Identifiability: nearest-(matched)-centroid classification of the
   withheld patient's cells using the LOPO-fitted centroids, compared
   against their `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R` full-dataset archetype label -- can the
   withheld patient's niches be recovered purely from everyone else's
   structure?
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_default_seed

SCALE_PREFIX = "radius_30.0um__"
RNG_SEED = get_default_seed()
N_INIT = 10


def match_centroids(lopo_centroids: np.ndarray, reference_centroids: np.ndarray) -> np.ndarray:
    """Pure, testable: optimal (Hungarian) matching of LOPO centroids to
    reference centroids by Euclidean distance. Returns `matching` such
    that `lopo_centroids[i]` corresponds to
    `reference_centroids[matching[i]]`."""
    cost = cdist(lopo_centroids, reference_centroids, metric="euclidean")
    row_idx, col_idx = linear_sum_assignment(cost)
    matching = np.empty(len(lopo_centroids), dtype=int)
    matching[row_idx] = col_idx
    return matching


def compute_centroid_cosine_similarity(
    lopo_centroids: np.ndarray, reference_centroids: np.ndarray, matching: np.ndarray
) -> np.ndarray:
    """Pure, testable: per-archetype cosine similarity between each LOPO
    centroid and its Hungarian-matched reference centroid."""
    matched_reference = reference_centroids[matching]
    dot = np.sum(lopo_centroids * matched_reference, axis=1)
    norm_lopo = np.linalg.norm(lopo_centroids, axis=1)
    norm_ref = np.linalg.norm(matched_reference, axis=1)
    return dot / (norm_lopo * norm_ref)


def assign_nearest_centroid(data: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    """Pure, testable: hard nearest-centroid (Euclidean) assignment."""
    dist = cdist(data, centroids, metric="euclidean")
    return np.argmin(dist, axis=1)


def build_lopo_stability(project_root: Path) -> dict:
    compositions_path = project_root / "data" / "derived" / "local_compositions.parquet"
    archetypes_path = project_root / "data" / "derived" / "neighbourhood_archetypes.parquet"
    centroids_path = project_root / "data" / "derived" / "neighbourhood_archetype_centroids.parquet"
    manifest_path = project_root / "metadata" / "sample_manifest.tsv"
    output_path = project_root / "data" / "derived" / "lopo_stability_results.parquet"

    for path, phase in [
        (
            compositions_path,
            "`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`",
        ),
        (
            archetypes_path,
            "`10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`",
        ),
        (
            centroids_path,
            "`10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`",
        ),
        (manifest_path, None),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found." + (f" Run {phase} first." if phase else ""))

    # `local_compositions.parquet` was written by Python (`10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py`):
    # pandas round-trips its named index (`cell_id`) as an index, not a
    # column -- unlike neighbourhood_archetypes.parquet, which was
    # written by R and materialises `cell_id` as a genuine column.
    compositions = pd.read_parquet(compositions_path)
    scale_cols = [c for c in compositions.columns if c.startswith(SCALE_PREFIX)]
    lineage_names = [c[len(SCALE_PREFIX) :] for c in scale_cols]

    archetypes = pd.read_parquet(archetypes_path).set_index("cell_id")
    reference_centroids_df = pd.read_parquet(centroids_path).set_index("archetype").sort_index()
    reference_centroids = reference_centroids_df[lineage_names].to_numpy()
    n_archetypes = len(reference_centroids_df)

    manifest = pd.read_csv(manifest_path, sep="\t")
    section_to_patient = manifest.set_index("section_id")["patient_id"]

    complete_mask = compositions[scale_cols].notna().all(axis=1)
    fitted_compositions = compositions.loc[complete_mask]
    fitted_patient = fitted_compositions["section_id"].map(section_to_patient)

    real_archetype = archetypes["archetype"].reindex(fitted_compositions.index).to_numpy()

    patients = sorted(fitted_patient.dropna().unique())
    fold_rows = []
    for patient in patients:
        held_out_mask = (fitted_patient == patient).to_numpy()
        train_mask = ~held_out_mask

        train_data = fitted_compositions.loc[train_mask, scale_cols].to_numpy()
        held_out_data = fitted_compositions.loc[held_out_mask, scale_cols].to_numpy()
        held_out_real_archetype = real_archetype[held_out_mask]

        km = KMeans(n_clusters=n_archetypes, n_init=N_INIT, random_state=RNG_SEED)
        km.fit(train_data)
        lopo_centroids = km.cluster_centers_

        matching = match_centroids(lopo_centroids, reference_centroids)
        cosine_similarity = compute_centroid_cosine_similarity(
            lopo_centroids, reference_centroids, matching
        )

        held_out_lopo_code = assign_nearest_centroid(held_out_data, lopo_centroids)
        # Map LOPO cluster codes (0-indexed) through the matching to
        # archetype ids (reference_centroids_df.index, e.g. 1..K).
        code_to_archetype_id = {
            code: reference_centroids_df.index[matched] for code, matched in enumerate(matching)
        }
        held_out_predicted_archetype = np.array(
            [code_to_archetype_id[c] for c in held_out_lopo_code]
        )

        identifiability_accuracy = float(
            np.mean(held_out_predicted_archetype == held_out_real_archetype)
        )

        for archetype_id, similarity in zip(
            reference_centroids_df.index[matching], cosine_similarity
        ):
            fold_rows.append(
                {
                    "patient_id": patient,
                    "n_held_out_cells": int(held_out_mask.sum()),
                    "n_train_cells": int(train_mask.sum()),
                    "archetype": archetype_id,
                    "centroid_cosine_similarity": float(similarity),
                    "identifiability_accuracy": identifiability_accuracy,
                }
            )

    result = pd.DataFrame(fold_rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    per_patient_accuracy = result.drop_duplicates("patient_id").set_index("patient_id")[
        "identifiability_accuracy"
    ]

    return {
        "n_patients": len(patients),
        "n_archetypes": n_archetypes,
        "mean_identifiability_accuracy": float(per_patient_accuracy.mean()),
        "min_identifiability_accuracy": float(per_patient_accuracy.min()),
        "mean_centroid_cosine_similarity": float(result["centroid_cosine_similarity"].mean()),
        "min_centroid_cosine_similarity": float(result["centroid_cosine_similarity"].min()),
        "output_path": str(output_path),
    }
