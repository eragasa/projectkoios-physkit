"""Verification of finite-difference Hamiltonian composition."""

import numpy as np

from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.discretization.cartesian_grids import UniformCartesianGrid1D
from projectkoios.physkit.discretization.intervals import DirichletInterval
from projectkoios.physkit.operators.operators_1d import (
    SecondOrderCentralDifferenceLaplacian1D,
)
from projectkoios.physkit.operators.qm import (
    FiniteDifferenceHamiltonian1D,
    SampledPotential1D,
    TiseKineticEnergy1D,
)
from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)


def test_adds_unitless_kinetic_and_potential() -> None:
    grid = UniformCartesianGrid1D(
        ScalarQuantity(-1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(0.5, Unitless()),
    )
    laplacian = SecondOrderCentralDifferenceLaplacian1D(
        DirichletInterval(
            grid,
            DirichletBoundaryCondition(ScalarQuantity(0.0, Unitless())),
        )
    )
    kinetic = TiseKineticEnergy1D(
        laplacian,
        ScalarQuantity(1.0, Unitless()),
        ScalarQuantity(1.0, Unitless()),
        Unitless(),
    )
    potential = SampledPotential1D(
        grid,
        VectorQuantity(np.array([0.125, 0.0, 0.125]), Unitless()),
    )

    matrix = FiniteDifferenceHamiltonian1D(kinetic, potential).matrix()

    np.testing.assert_array_equal(
        matrix.to_dense().magnitude,
        [[4.125, -2.0, 0.0], [-2.0, 4.0, -2.0], [0.0, -2.0, 4.125]],
    )


def test_converts_potential_to_selected_energy_unit() -> None:
    grid = UniformCartesianGrid1D(
        ScalarQuantity(0.0, PhysicalUnit("angstrom")),
        ScalarQuantity(10.0, PhysicalUnit("angstrom")),
        ScalarQuantity(2.5, PhysicalUnit("angstrom")),
    )
    laplacian = SecondOrderCentralDifferenceLaplacian1D(
        DirichletInterval(
            grid,
            DirichletBoundaryCondition(
                ScalarQuantity(0.0, PhysicalUnit("angstrom ** -0.5"))
            ),
        )
    )
    kinetic = TiseKineticEnergy1D(
        laplacian,
        ScalarQuantity(
            6.582_119_569e-4,
            PhysicalUnit("electron_volt * picosecond"),
        ),
        ScalarQuantity(5.485_799_090_65e-4, PhysicalUnit("dalton")),
        PhysicalUnit("electron_volt"),
    )
    potential = SampledPotential1D(
        grid,
        VectorQuantity(
            np.array([100.0, 200.0, 300.0]),
            PhysicalUnit("millielectron_volt"),
        ),
    )

    matrix = FiniteDifferenceHamiltonian1D(kinetic, potential).matrix()
    expected = MODEL_SYSTEM_UNIT_CONVERTER.convert_vector(
        potential.values,
        PhysicalUnit("electron_volt"),
    ).magnitude

    np.testing.assert_allclose(
        np.diag(matrix.to_dense().magnitude)
        - np.diag(kinetic.matrix().to_dense().magnitude),
        expected,
        rtol=2.0e-15,
        atol=0.0,
    )
    assert matrix.unit == PhysicalUnit("electron_volt")
