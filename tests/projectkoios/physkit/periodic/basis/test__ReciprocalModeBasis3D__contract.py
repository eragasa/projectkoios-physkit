"""Contract tests for the canonical three-dimensional reciprocal-mode basis."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis3D
from projectkoios.physkit.periodic.lattice import (
    ReciprocalLattice2D,
    ReciprocalLattice3D,
)


def _reciprocal_lattice() -> ReciprocalLattice3D:
    return ReciprocalLattice3D(
        b1=np.asarray([2.0, 0.0, 0.0]),
        b2=np.asarray([0.5, 3.0, 0.0]),
        b3=np.asarray([-0.25, 0.75, 4.0]),
    )


def test__init__retains_explicit_mode_order_and_constructs_vectors() -> None:
    """Explicit mode ordering survives canonical lattice-vector construction."""
    reciprocal = _reciprocal_lattice()
    source_modes = np.asarray([(1, -1, 0), (0, 0, 0), (-2, 1, 3)])

    basis = ReciprocalModeBasis3D(reciprocal, source_modes)

    assert basis.reciprocal_lattice is reciprocal
    assert basis.size == 3
    np.testing.assert_array_equal(basis.modes, source_modes)
    np.testing.assert_allclose(
        basis.vectors,
        [[1.5, -3.0, 0.0], [0.0, 0.0, 0.0], [-4.25, 5.25, 12.0]],
    )
    assert basis.modes.flags.writeable is False
    assert basis.vectors.flags.writeable is False
    assert source_modes.flags.writeable is True


def test__from_mode_limit__uses_first_outer_and_third_inner_order() -> None:
    """The cubic factory exposes one deterministic tensor-product ordering."""
    basis = ReciprocalModeBasis3D.from_mode_limit(_reciprocal_lattice(), 1)

    assert basis.size == 27
    np.testing.assert_array_equal(
        basis.modes[:4],
        [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1)],
    )
    np.testing.assert_array_equal(basis.modes[-1], (1, 1, 1))


def test__from_mode_limit__supports_the_single_constant_mode() -> None:
    """Cutoff zero represents exactly the reciprocal origin."""
    basis = ReciprocalModeBasis3D.from_mode_limit(_reciprocal_lattice(), 0)

    assert basis.size == 1
    np.testing.assert_array_equal(basis.modes, [(0, 0, 0)])
    np.testing.assert_array_equal(basis.vectors, np.zeros((1, 3)))


def test__shifted_wavevectors__adds_one_cartesian_vector() -> None:
    """A finite Cartesian carrier shifts every represented reciprocal vector."""
    basis = ReciprocalModeBasis3D(
        _reciprocal_lattice(),
        [(0, 0, 0), (1, -1, 2)],
    )
    shift = np.asarray([0.25, -0.5, 0.75])

    np.testing.assert_allclose(basis.shifted_wavevectors(shift), basis.vectors + shift)


def test__init__rejects_noncanonical_lattice_and_invalid_inventory() -> None:
    """Nominal ownership, dimensional shape, integer type, and uniqueness are strict."""
    reciprocal = _reciprocal_lattice()
    reciprocal_2d = ReciprocalLattice2D(
        b1=np.asarray([1.0, 0.0]),
        b2=np.asarray([0.0, 1.0]),
    )

    with pytest.raises(TypeError, match="must be ReciprocalLattice3D"):
        ReciprocalModeBasis3D(reciprocal_2d, [(0, 0, 0)])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="cannot be empty"):
        ReciprocalModeBasis3D(reciprocal, [])
    with pytest.raises(ValueError, match="must have shape"):
        ReciprocalModeBasis3D(reciprocal, [(0, 0)])  # type: ignore[list-item]
    with pytest.raises(TypeError, match="must contain integers"):
        ReciprocalModeBasis3D(reciprocal, [(0.0, 0.0, 0.0)])  # type: ignore[list-item]
    with pytest.raises(ValueError, match="must be unique"):
        ReciprocalModeBasis3D(reciprocal, [(0, 0, 0), (0, 0, 0)])
    with pytest.raises(ValueError, match="signed 64-bit"):
        ReciprocalModeBasis3D(
            reciprocal,
            np.asarray([[np.iinfo(np.uint64).max, 0, 0]], dtype=np.uint64),
        )


def test__from_mode_limit__and_shift_reject_invalid_values() -> None:
    """Cutoff and shift validation reject ambiguous or nonfinite coordinates."""
    reciprocal = _reciprocal_lattice()

    with pytest.raises(TypeError, match="built-in integer"):
        ReciprocalModeBasis3D.from_mode_limit(reciprocal, True)
    with pytest.raises(ValueError, match="must be nonnegative"):
        ReciprocalModeBasis3D.from_mode_limit(reciprocal, -1)

    basis = ReciprocalModeBasis3D.from_mode_limit(reciprocal, 0)
    with pytest.raises(ValueError, match="shape"):
        basis.shifted_wavevectors([0.0, 0.0])
    with pytest.raises(TypeError, match="real numeric"):
        basis.shifted_wavevectors([0.0, 0.0, 1.0j])
    with pytest.raises(ValueError, match="finite"):
        basis.shifted_wavevectors([0.0, np.inf, 0.0])
