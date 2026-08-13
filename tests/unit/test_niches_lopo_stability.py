"""Unit tests for xenium_tcr_ecology.niches.lopo_stability (`10_niche_and_ecosystem_discovery/07_leave_one_patient_out_niche_stability.py`)."""

from __future__ import annotations

import numpy as np

from xenium_tcr_ecology.niches.lopo_stability import (
    assign_nearest_centroid,
    compute_centroid_cosine_similarity,
    match_centroids,
)


class TestMatchCentroids:
    def test_identity_when_already_aligned(self):
        reference = np.array([[0.0, 0.0], [10.0, 10.0], [20.0, 0.0]])
        lopo = np.array([[0.1, 0.1], [10.1, 9.9], [19.9, 0.1]])
        matching = match_centroids(lopo, reference)
        assert list(matching) == [0, 1, 2]

    def test_recovers_permutation(self):
        reference = np.array([[0.0, 0.0], [10.0, 10.0], [20.0, 0.0]])
        # LOPO centroids in a shuffled order relative to reference.
        lopo = np.array([[20.1, 0.1], [0.1, -0.1], [9.9, 10.1]])
        matching = match_centroids(lopo, reference)
        # lopo[0] (near ref 2), lopo[1] (near ref 0), lopo[2] (near ref 1)
        assert list(matching) == [2, 0, 1]


class TestComputeCentroidCosineSimilarity:
    def test_identical_matched_centroids_give_similarity_one(self):
        reference = np.array([[1.0, 0.0], [0.0, 1.0]])
        lopo = np.array([[2.0, 0.0], [0.0, 3.0]])  # same direction, different magnitude
        matching = np.array([0, 1])
        result = compute_centroid_cosine_similarity(lopo, reference, matching)
        assert np.allclose(result, [1.0, 1.0])

    def test_orthogonal_matched_centroids_give_similarity_zero(self):
        reference = np.array([[1.0, 0.0]])
        lopo = np.array([[0.0, 1.0]])
        matching = np.array([0])
        result = compute_centroid_cosine_similarity(lopo, reference, matching)
        assert np.allclose(result, [0.0])


class TestAssignNearestCentroid:
    def test_assigns_to_the_closest_centroid(self):
        centroids = np.array([[0.0, 0.0], [10.0, 10.0]])
        data = np.array([[0.5, 0.5], [9.5, 9.5], [0.1, -0.1]])
        result = assign_nearest_centroid(data, centroids)
        assert list(result) == [0, 1, 0]
