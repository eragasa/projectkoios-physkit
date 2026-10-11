from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D


def test_stores_primitive_vectors_as_immutable_columns() -> None:
    a1 = np.array([2.0, 1.0])
    a2 = np.array([1.0, 3.0])

    lattice = DirectLattice2D(a1=a1, a2=a2)
    a1[0] = 99.0

    np.testing.assert_array_equal(
        lattice.A,
        np.array([[2.0, 1.0], [1.0, 3.0]]),
    )
    assert lattice.A.flags.writeable is False


def test_rejects_linearly_dependent_vectors() -> None:
    with pytest.raises(ValueError, match="linearly independent"):
        DirectLattice2D(
            a1=np.array([1.0, 0.0]),
            a2=np.array([2.0, 0.0]),
        )


def test__init__uses_scale_invariant_independence_check() -> None:
    lattice = DirectLattice2D(
        a1=np.array([3.0e-10, 0.0]),
        a2=np.array([1.0e-10, 2.0e-10]),
    )

    np.testing.assert_array_equal(
        lattice.A,
        np.array([[3.0e-10, 1.0e-10], [0.0, 2.0e-10]]),
    )
