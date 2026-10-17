"""Contract tests for the species-neutral diamond-decorated Bain path."""

from __future__ import annotations

import numpy as np
import pytest

from projectkoios.physkit.periodic.lattice import (
    DiamondDecoratedBainPath3D,
    DiamondDecoratedBainPathPoint3D,
    UnitVolumeBainPath3D,
)
from projectkoios.physkit.units import MatrixQuantity, Unitless


def test__evaluate__retains_bain_and_fcc_fractional_basis_conventions() -> None:
    """The decorated path owns both exact fractional-basis representations."""
    point = DiamondDecoratedBainPath3D().evaluate(0.0)

    assert type(point) is DiamondDecoratedBainPathPoint3D
    assert type(point.fractional_basis_bain) is MatrixQuantity
    assert type(point.fractional_basis_bain.unit) is Unitless
    assert type(point.fractional_basis_fcc.unit) is Unitless
    np.testing.assert_array_equal(
        point.fractional_basis_bain.magnitude,
        [[0.0, 0.0, 0.0], [0.5, -0.25, 0.25]],
    )
    np.testing.assert_array_equal(
        point.fractional_basis_fcc.magnitude,
        [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    )


def test__evaluate__affinely_carries_basis_through_bain_lattice() -> None:
    """Cartesian decorations derive from the current lattice and fixed fractions."""
    path = DiamondDecoratedBainPath3D()

    for parameter in (0.0, 0.25, 0.5, 0.75, 1.0):
        point = path.evaluate(parameter)
        expected = point.fractional_basis_bain.magnitude @ point.direct_lattice.A.T
        np.testing.assert_array_equal(point.cartesian_basis.magnitude, expected)
        assert type(point.cartesian_basis.unit) is Unitless
        assert point.primitive_volume == pytest.approx(1.0, rel=0.0, abs=5.0e-16)


def test__point__owns_correlated_bravais_path_state() -> None:
    """The decorated record copies the complete Bravais path point."""
    source = UnitVolumeBainPath3D().evaluate(0.5)

    decorated = DiamondDecoratedBainPathPoint3D(source)

    assert decorated.bain_point is not source
    assert decorated.direct_lattice is decorated.bain_point.direct_lattice
    assert decorated.direct_lattice is not source.direct_lattice
    np.testing.assert_array_equal(decorated.direct_lattice.A, source.direct_lattice.A)


def test__point__rejects_non_bain_state() -> None:
    """Diamond decoration requires the exact correlated Bravais owner."""
    with pytest.raises(TypeError, match="must be UnitVolumeBainPathPoint3D"):
        DiamondDecoratedBainPathPoint3D(object())  # type: ignore[arg-type]
