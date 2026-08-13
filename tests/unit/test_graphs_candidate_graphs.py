"""Unit tests for xenium_tcr_ecology.graphs.candidate_graphs (`09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`)."""

from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon

from xenium_tcr_ecology.graphs.candidate_graphs import (
    build_boundary_contact_graph,
    build_delaunay_graph,
    build_knn_graph,
    build_radius_graph,
)


def _square(x0, y0, side=1.0):
    return Polygon([(x0, y0), (x0 + side, y0), (x0 + side, y0 + side), (x0, y0 + side)])


class TestBuildBoundaryContactGraph:
    def test_touching_polygons_are_connected(self):
        polys = np.array([_square(0, 0), _square(1, 0)], dtype=object)  # share an edge
        graph = build_boundary_contact_graph(polys)
        assert graph[0, 1] != 0
        assert graph[1, 0] != 0

    def test_non_touching_polygons_are_not_connected(self):
        polys = np.array([_square(0, 0), _square(100, 100)], dtype=object)
        graph = build_boundary_contact_graph(polys)
        assert graph[0, 1] == 0

    def test_no_self_loops(self):
        polys = np.array([_square(0, 0), _square(1, 0)], dtype=object)
        graph = build_boundary_contact_graph(polys)
        assert graph[0, 0] == 0
        assert graph[1, 1] == 0

    def test_empty_input_returns_empty_graph(self):
        graph = build_boundary_contact_graph(np.array([], dtype=object))
        assert graph.shape == (0, 0)

    def test_small_real_world_gap_is_connected_within_tolerance(self):
        # Regression test: real Xenium cell-boundary polygons are not
        # seamlessly touching -- most adjacent cells have a small nonzero
        # gap (median ~0.08um, checked directly on real data). A strict
        # `intersects`-only definition gave an implausible mean degree of
        # 0.18 on a real 25,964-cell section (vs. ~6 for Delaunay on the
        # same cells).
        polys = np.array([_square(0, 0), _square(1.3, 0)], dtype=object)  # 0.3um gap
        graph = build_boundary_contact_graph(polys, tolerance_um=0.5)
        assert graph[0, 1] != 0

    def test_gap_beyond_tolerance_is_not_connected(self):
        polys = np.array([_square(0, 0), _square(5.0, 0)], dtype=object)  # 4.0um gap
        graph = build_boundary_contact_graph(polys, tolerance_um=0.5)
        assert graph[0, 1] == 0


class TestBuildRadiusGraph:
    def test_cells_within_radius_are_connected(self):
        coords = np.array([[0.0, 0.0], [5.0, 0.0], [1000.0, 0.0]])
        graph = build_radius_graph(coords, radius_um=10.0)
        assert graph[0, 1] != 0
        assert graph[0, 2] == 0


class TestBuildKnnGraph:
    def test_connects_to_nearest_neighbors(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1000.0, 0.0]])
        graph = build_knn_graph(coords, n_neighs=2)
        assert graph[0, 1] != 0

    def test_handles_fewer_cells_than_requested_neighbors(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0]])
        graph = build_knn_graph(coords, n_neighs=10)
        assert graph.shape == (2, 2)

    def test_single_cell_returns_empty_graph(self):
        coords = np.array([[0.0, 0.0]])
        graph = build_knn_graph(coords, n_neighs=6)
        assert graph.nnz == 0


class TestBuildDelaunayGraph:
    def test_fewer_than_four_points_returns_empty_graph(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        graph = build_delaunay_graph(coords)
        assert graph.nnz == 0
        assert graph.shape == (3, 3)

    def test_four_points_produces_a_real_triangulation(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        graph = build_delaunay_graph(coords)
        assert graph.nnz > 0
