"""Contract tests for ambient three-dimensional 2D direct lattices."""

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice.embedding import EmbeddedDirectLattice2D
from projectkoios.physkit.units import MatrixQuantity, PhysicalUnit, Unitless


def test__init__retains_immutable_ambient_orientation_and_unit() -> None:
    values = np.array(
        [
            [2.0, 0.4],
            [0.3, 1.7],
            [1.1, -0.2],
        ]
    )
    embedded = EmbeddedDirectLattice2D(
        ambient_basis=MatrixQuantity(values, PhysicalUnit("angstrom"))
    )
    values[0, 0] = 99.0

    np.testing.assert_array_equal(
        embedded.ambient_basis.magnitude,
        np.array(
            [
                [2.0, 0.4],
                [0.3, 1.7],
                [1.1, -0.2],
            ]
        ),
    )
    assert embedded.ambient_basis.magnitude.flags.writeable is False
    assert embedded.basis_unit == PhysicalUnit("angstrom")


def test__init__uses_scale_invariant_full_column_rank() -> None:
    direction_basis = np.array(
        [
            [1.0, 0.2],
            [0.1, 1.0],
            [0.3, -0.4],
        ]
    )

    small = EmbeddedDirectLattice2D(
        ambient_basis=MatrixQuantity(
            1.0e-10 * direction_basis,
            PhysicalUnit("meter"),
        )
    )
    independently_scaled = EmbeddedDirectLattice2D(
        ambient_basis=MatrixQuantity(
            direction_basis * np.array([1.0e-7, 1.0e7]),
            Unitless(),
        )
    )

    assert small.ambient_basis.magnitude.shape == (3, 2)
    assert independently_scaled.ambient_basis.magnitude.shape == (3, 2)


def test__init__rejects_wrong_shape_rank_and_unit() -> None:
    with pytest.raises(ValueError, match=r"shape \(3, 2\)"):
        EmbeddedDirectLattice2D(
            ambient_basis=MatrixQuantity(np.identity(2), Unitless())
        )

    with pytest.raises(ValueError, match="full column rank"):
        EmbeddedDirectLattice2D(
            ambient_basis=MatrixQuantity(
                np.array(
                    [
                        [1.0, 2.0],
                        [2.0, 4.0],
                        [3.0, 6.0],
                    ]
                ),
                Unitless(),
            )
        )

    with pytest.raises(ValueError, match="full column rank"):
        EmbeddedDirectLattice2D(
            ambient_basis=MatrixQuantity(
                np.array(
                    [
                        [1.0, 1.0],
                        [0.0, 1.0e-10],
                        [0.0, 0.0],
                    ]
                ),
                Unitless(),
            )
        )

    with pytest.raises(ValueError, match="compatible with length"):
        EmbeddedDirectLattice2D(
            ambient_basis=MatrixQuantity(
                np.array(
                    [
                        [1.0, 0.0],
                        [0.0, 1.0],
                        [0.0, 0.0],
                    ]
                ),
                PhysicalUnit("electron_volt"),
            )
        )
