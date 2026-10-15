"""Verification of the represented one-dimensional Laplacian."""

import numpy as np
import pytest

from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.discretization.cartesian_grids import UniformCartesianGrid1D
from projectkoios.physkit.discretization.intervals import DirichletInterval
from projectkoios.physkit.operators.operators_1d import (
    SecondOrderCentralDifferenceLaplacian1D,
    UnitAwareMatrixOperator,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, Unitless


def test_returns_centered_unitless_dirichlet_stencil() -> None:
    assert issubclass(
        SecondOrderCentralDifferenceLaplacian1D,
        UnitAwareMatrixOperator,
    )
    grid = UniformCartesianGrid1D(
        ScalarQuantity(-1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(0.5, Unitless()),
    )
    matrix = SecondOrderCentralDifferenceLaplacian1D(
        DirichletInterval(
            grid,
            DirichletBoundaryCondition(ScalarQuantity(0.0, Unitless())),
        )
    ).matrix()

    np.testing.assert_array_equal(
        matrix.to_dense().magnitude,
        [[-8.0, 4.0, 0.0], [4.0, -8.0, 4.0], [0.0, 4.0, -8.0]],
    )
    assert isinstance(matrix.unit, Unitless)
    assert not matrix.data.flags.writeable


def test_preserves_selected_physical_length_scale() -> None:
    grid = UniformCartesianGrid1D(
        ScalarQuantity(0.0, PhysicalUnit("angstrom")),
        ScalarQuantity(10.0, PhysicalUnit("angstrom")),
        ScalarQuantity(2.5, PhysicalUnit("angstrom")),
    )
    matrix = SecondOrderCentralDifferenceLaplacian1D(
        DirichletInterval(
            grid,
            DirichletBoundaryCondition(
                ScalarQuantity(0.0, PhysicalUnit("angstrom ** -0.5"))
            ),
        )
    ).matrix()

    np.testing.assert_allclose(
        matrix.to_dense().magnitude,
        np.array([[-2.0, 1.0, 0.0], [1.0, -2.0, 1.0], [0.0, 1.0, -2.0]]) / 2.5**2,
        rtol=2.0e-16,
        atol=0.0,
    )
    assert matrix.unit == PhysicalUnit("1 / (angstrom) ** 2")


def test_rejects_nonhomogeneous_boundary() -> None:
    grid = UniformCartesianGrid1D(
        ScalarQuantity(-1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(0.5, Unitless()),
    )
    with pytest.raises(ValueError, match="requires homogeneous"):
        SecondOrderCentralDifferenceLaplacian1D(
            DirichletInterval(
                grid,
                DirichletBoundaryCondition(ScalarQuantity(1.0, Unitless())),
            )
        )
