"""Verification of represented Schrödinger kinetic energy."""

import numpy as np

from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.discretization.cartesian_grids import UniformCartesianGrid1D
from projectkoios.physkit.discretization.intervals import DirichletInterval
from projectkoios.physkit.operators.operators_1d import (
    SecondOrderCentralDifferenceLaplacian1D,
)
from projectkoios.physkit.operators.qm import TiseKineticEnergy1D
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, Unitless


def laplacian(
    length: float,
    spacing: float,
    unit: PhysicalUnit | Unitless,
) -> SecondOrderCentralDifferenceLaplacian1D:
    boundary_unit = (
        Unitless()
        if isinstance(unit, Unitless)
        else PhysicalUnit(f"({unit.expression}) ** -0.5")
    )
    return SecondOrderCentralDifferenceLaplacian1D(
        DirichletInterval(
            UniformCartesianGrid1D(
                ScalarQuantity(0.0, unit),
                ScalarQuantity(length, unit),
                ScalarQuantity(spacing, unit),
            ),
            DirichletBoundaryCondition(ScalarQuantity(0.0, boundary_unit)),
        )
    )


def test_applies_negative_half_unitless_laplacian_scale() -> None:
    matrix = TiseKineticEnergy1D(
        laplacian(2.0, 0.5, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        Unitless(),
    ).matrix()

    np.testing.assert_array_equal(
        matrix.to_dense().magnitude,
        [[4.0, -2.0, 0.0], [-2.0, 4.0, -2.0], [0.0, -2.0, 4.0]],
    )


def test_returns_well_scaled_metal_energy_matrix() -> None:
    matrix = TiseKineticEnergy1D(
        laplacian(10.0, 2.5, PhysicalUnit("angstrom")),
        ScalarQuantity(
            6.582_119_569e-4,
            PhysicalUnit("electron_volt * picosecond"),
        ),
        ScalarQuantity(
            5.485_799_090_65e-4,
            PhysicalUnit("dalton"),
        ),
        PhysicalUnit("electron_volt"),
    ).matrix()

    assert matrix.unit == PhysicalUnit("electron_volt")
    assert np.all(np.diag(matrix.to_dense().magnitude) > 0.0)
    assert np.all(np.diag(matrix.to_dense().magnitude, 1) < 0.0)
    assert np.max(np.abs(matrix.data)) < 10.0
