"""Contract tests for the exact FCC--Bain coordinate bridge."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice import FccBainCoordinateBridge3D


def test__matrices__are_exact_unimodular_inverses() -> None:
    """The declared integer maps compose to identity in both directions."""
    bridge = FccBainCoordinateBridge3D()

    np.testing.assert_array_equal(
        bridge.bain_to_fcc_matrix,
        [[1, 1, 0], [0, 0, 1], [0, -1, 0]],
    )
    np.testing.assert_array_equal(
        bridge.fcc_to_bain_matrix,
        [[1, 0, 1], [0, 0, -1], [0, 1, 0]],
    )
    np.testing.assert_array_equal(
        bridge.bain_to_fcc_matrix @ bridge.fcc_to_bain_matrix,
        np.eye(3, dtype=np.int64),
    )
    assert bridge.determinant == 1
    assert bridge.bain_to_fcc_matrix.flags.writeable is False
    assert bridge.fcc_to_bain_matrix.flags.writeable is False


def test__fractional_maps__recover_diamond_offset_and_round_trip() -> None:
    """The nonzero Bain decoration maps to the conventional FCC diamond offset."""
    bridge = FccBainCoordinateBridge3D()
    bain = np.asarray([[0.0, 0.0, 0.0], [0.5, -0.25, 0.25]])

    fcc = bridge.bain_to_fcc_fractional(bain)

    np.testing.assert_array_equal(fcc, [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
    np.testing.assert_array_equal(bridge.fcc_to_bain_fractional(fcc), bain)


def test__integer_maps__round_trip_and_remain_bijective_modulo_finite_sizes() -> None:
    """Unimodularity preserves every finite periodic quotient inventory."""
    bridge = FccBainCoordinateBridge3D()
    indices = np.asarray([(2, -1, 3), (0, 0, 0), (-4, 5, 1)], dtype=np.int64)

    mapped = bridge.bain_to_fcc_indices(indices)

    np.testing.assert_array_equal(bridge.fcc_to_bain_indices(mapped), indices)
    for size in range(2, 9):
        images = {
            tuple(bridge.bain_to_fcc_indices(np.asarray(index)) % size)
            for index in itertools.product(range(size), repeat=3)
        }
        assert len(images) == size**3


def test__maps__reject_wrong_shape_type_range_or_nonfinite_coordinates() -> None:
    """Coordinate kind and final dimension remain explicit at the public boundary."""
    bridge = FccBainCoordinateBridge3D()

    with pytest.raises(ValueError, match="final dimension 3"):
        bridge.bain_to_fcc_indices(np.asarray([1, 2]))
    with pytest.raises(TypeError, match="must contain integers"):
        bridge.bain_to_fcc_indices(np.asarray([1.0, 2.0, 3.0]))
    with pytest.raises(ValueError, match="signed 64-bit"):
        bridge.bain_to_fcc_indices(
            np.asarray([np.iinfo(np.uint64).max, 0, 0], dtype=np.uint64)
        )
    with pytest.raises(TypeError, match="real numeric"):
        bridge.bain_to_fcc_fractional(np.asarray([0.0, 0.0, 1.0j]))
    with pytest.raises(ValueError, match="finite"):
        bridge.bain_to_fcc_fractional(np.asarray([0.0, np.inf, 0.0]))
