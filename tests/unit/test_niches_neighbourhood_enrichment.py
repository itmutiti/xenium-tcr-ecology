"""Unit tests for xenium_tcr_ecology.niches.neighbourhood_enrichment (`10_niche_and_ecosystem_discovery/00_compute_cell_type_neighbourhood_enrichment.py`)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from xenium_tcr_ecology.niches.neighbourhood_enrichment import compute_section_enrichment


def _make_two_cluster_data(rng, n_per_cluster=60):
    # Cluster A cells connect densely to each other and rarely to cluster
    # B -- a real, self-avoiding neighbourhood structure to test against.
    n = n_per_cluster * 2
    lineage = ["A"] * n_per_cluster + ["B"] * n_per_cluster
    node_metadata = pd.DataFrame({"final_lineage": lineage}, index=[f"c{i}" for i in range(n)])

    rows, cols = [], []
    for i in range(n_per_cluster):
        for j in range(i + 1, n_per_cluster):
            if rng.random() < 0.3:
                rows += [i, j]
                cols += [j, i]
    for i in range(n_per_cluster, n):
        for j in range(i + 1, n):
            if rng.random() < 0.3:
                rows += [i, j]
                cols += [j, i]
    # A handful of real cross-cluster edges.
    for _ in range(5):
        a, b = rng.integers(0, n_per_cluster), rng.integers(n_per_cluster, n)
        rows += [a, b]
        cols += [b, a]

    data = [1.0] * len(rows)
    graph = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
    return node_metadata, graph


class TestComputeSectionEnrichment:
    def test_self_avoiding_clusters_show_depletion_between_groups(self):
        rng = np.random.default_rng(0)
        node_metadata, graph = _make_two_cluster_data(rng)
        result = compute_section_enrichment(node_metadata, graph, n_perms=99, seed=0)
        cross = result[
            ((result["lineage_a"] == "A") & (result["lineage_b"] == "B"))
            | ((result["lineage_a"] == "B") & (result["lineage_b"] == "A"))
        ]
        same = result[(result["lineage_a"] == "A") & (result["lineage_b"] == "A")]
        assert cross.iloc[0]["zscore"] < same.iloc[0]["zscore"]

    def test_returns_empty_when_fewer_than_two_eligible_lineages(self):
        node_metadata = pd.DataFrame(
            {"final_lineage": ["A"] * 5 + ["B"] * 2}, index=[f"c{i}" for i in range(7)]
        )
        graph = sparse.csr_matrix((7, 7))
        result = compute_section_enrichment(node_metadata, graph, n_perms=99)
        assert len(result) == 0

    def test_output_has_no_duplicate_unordered_pairs(self):
        rng = np.random.default_rng(1)
        node_metadata, graph = _make_two_cluster_data(rng)
        result = compute_section_enrichment(node_metadata, graph, n_perms=99, seed=1)
        pairs = set(zip(result["lineage_a"], result["lineage_b"]))
        assert ("A", "B") not in pairs or ("B", "A") not in pairs
