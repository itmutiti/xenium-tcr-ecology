"""Unit tests for xenium_tcr_ecology.niches.ecosystem_metrics (`10_niche_and_ecosystem_discovery/05_quantify_ecosystem_abundance_and_topology.py`)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from xenium_tcr_ecology.niches.ecosystem_metrics import (
    compute_effective_n_patches,
    compute_interface_edge_count,
    compute_isoperimetric_quotient,
    compute_mixing_index,
)


class TestComputeEffectiveNPatches:
    def test_one_dominant_patch_plus_fragments_has_low_effective_count(self):
        # One big patch (90 cells) + 9 singleton fragments -- much closer
        # to 1 "effective" patch than to 10 equally-sized ones.
        sizes = np.array([90] + [1] * 9)
        result = compute_effective_n_patches(sizes)
        assert 1.0 < result < 2.0

    def test_equally_sized_patches_give_effective_count_equal_to_n(self):
        sizes = np.array([10, 10, 10, 10])
        result = compute_effective_n_patches(sizes)
        assert result == 4.0

    def test_single_patch_gives_one(self):
        assert compute_effective_n_patches(np.array([42])) == 1.0


class TestComputeInterfaceEdgeCount:
    def test_no_crossing_edges_when_all_same_code(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        codes = np.array([0, 0])
        assert compute_interface_edge_count(graph, codes, this_code=0) == 0

    def test_single_crossing_edge_counted_once(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        codes = np.array([0, 1])
        assert compute_interface_edge_count(graph, codes, this_code=0) == 1
        assert compute_interface_edge_count(graph, codes, this_code=1) == 1

    def test_interior_edges_of_a_third_code_are_not_counted(self):
        # Path 0-1-2, codes [0, 1, 2] -- edge (0,1) crosses for code 0 and
        # code 1; edge (1,2) crosses for code 1 and code 2; code 0's count
        # must not include the (1,2) edge it has no endpoint in.
        rows = [0, 1]
        cols = [1, 2]
        graph = sparse.csr_matrix(([1.0] * 4, (rows + cols, cols + rows)), shape=(3, 3))
        codes = np.array([0, 1, 2])
        assert compute_interface_edge_count(graph, codes, this_code=0) == 1


class TestComputeMixingIndex:
    def test_fully_segregated_block_has_zero_mixing(self):
        # Two disjoint same-code pairs -- every neighbour is same code.
        graph = sparse.csr_matrix(
            ([1.0, 1.0, 1.0, 1.0], ([0, 1, 2, 3], [1, 0, 3, 2])), shape=(4, 4)
        )
        codes = np.array([0, 0, 0, 0])
        assert compute_mixing_index(graph, codes, this_code=0) == 0.0

    def test_fully_interspersed_pair_has_full_mixing(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        codes = np.array([0, 1])
        assert compute_mixing_index(graph, codes, this_code=0) == 1.0

    def test_no_cells_of_this_code_gives_nan(self):
        graph = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))
        codes = np.array([1, 1])
        assert np.isnan(compute_mixing_index(graph, codes, this_code=0))


class TestComputeIsoperimetricQuotient:
    def test_square_is_less_compact_than_circle_but_well_defined(self):
        # Unit square corners.
        x = np.array([0.0, 1.0, 1.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        iq = compute_isoperimetric_quotient(x, y)
        assert iq is not None
        assert 0.7 < iq < 0.9  # a square's true IQ is pi/4 ~= 0.785

    def test_fewer_than_three_points_returns_none(self):
        assert compute_isoperimetric_quotient(np.array([0.0, 1.0]), np.array([0.0, 1.0])) is None

    def test_collinear_points_return_none(self):
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.0, 0.0, 0.0])
        assert compute_isoperimetric_quotient(x, y) is None

    def test_near_circular_hexagon_has_high_iq(self):
        angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        x = np.cos(angles)
        y = np.sin(angles)
        iq = compute_isoperimetric_quotient(x, y)
        assert iq is not None
        assert iq > 0.9
