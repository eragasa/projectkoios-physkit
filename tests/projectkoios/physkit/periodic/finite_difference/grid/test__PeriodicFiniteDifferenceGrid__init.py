"""Contract tests for endpoint-excluded multidimensional periodic grids."""

import numpy as np
import pytest

from projectkoios.physkit.periodic.finite_difference.grid import (
    PeriodicFiniteDifferenceGrid,
)
from projectkoios.physkit.periodic.lattice.finite_domain import (
    FinitePeriodicDomain,
    LatticeDimension,
)
from projectkoios.physkit.periodic.lattice.lattice2d import DirectLattice2D
from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.periodic.lattice.metric import DirectLatticeMetric
from projectkoios.physkit.units import PhysicalUnit, Unitless


def _metric_2d() -> DirectLatticeMetric:
    return DirectLatticeMetric(
        DirectLattice2D(
            a1=np.array([2.0, 0.0]),
            a2=np.array([0.5, 1.5]),
        )
    )


def test__init__retains_bounded_endpoint_excluded_grid_context() -> None:
    grid = PeriodicFiniteDifferenceGrid(
        metric=_metric_2d(),
        domain=FinitePeriodicDomain(
            LatticeDimension.TWO,
            (4, 5),
        ),
        direct_basis_unit=Unitless(),
    )

    assert grid.dimension is LatticeDimension.TWO
    assert grid.shape == (4, 5)
    assert grid.point_count == 20
    assert grid.fractional_spacings == (0.25, 0.2)
    assert isinstance(grid.laplacian_unit, Unitless)


def test__init__accepts_three_dimensional_physical_length_basis() -> None:
    metric = DirectLatticeMetric(
        DirectLattice3D(
            a1=np.array([3.0, 0.0, 0.0]),
            a2=np.array([0.5, 2.0, 0.0]),
            a3=np.array([0.25, 0.4, 1.5]),
        )
    )
    grid = PeriodicFiniteDifferenceGrid(
        metric=metric,
        domain=FinitePeriodicDomain(
            LatticeDimension.THREE,
            (3, 4, 5),
        ),
        direct_basis_unit=PhysicalUnit("angstrom"),
    )

    assert grid.fractional_spacings == (1.0 / 3.0, 0.25, 0.2)
    assert grid.laplacian_unit == PhysicalUnit("1 / (angstrom) ** 2")


def test__init__rejects_dimension_mismatch_and_one_dimensional_domain() -> None:
    with pytest.raises(ValueError, match="dimensions must agree"):
        PeriodicFiniteDifferenceGrid(
            metric=_metric_2d(),
            domain=FinitePeriodicDomain(
                LatticeDimension.THREE,
                (3, 3, 3),
            ),
            direct_basis_unit=Unitless(),
        )

    with pytest.raises(ValueError, match="dimension must be two or three"):
        PeriodicFiniteDifferenceGrid(
            metric=_metric_2d(),
            domain=FinitePeriodicDomain(
                LatticeDimension.ONE,
                (3,),
            ),
            direct_basis_unit=Unitless(),
        )


def test__init__rejects_extent_below_centered_stencil_minimum() -> None:
    with pytest.raises(ValueError, match="at least three"):
        PeriodicFiniteDifferenceGrid(
            metric=_metric_2d(),
            domain=FinitePeriodicDomain(
                LatticeDimension.TWO,
                (2, 4),
            ),
            direct_basis_unit=Unitless(),
        )


def test__init__rejects_physical_nonlength_basis_unit() -> None:
    with pytest.raises(ValueError, match="compatible with length"):
        PeriodicFiniteDifferenceGrid(
            metric=_metric_2d(),
            domain=FinitePeriodicDomain(
                LatticeDimension.TWO,
                (3, 3),
            ),
            direct_basis_unit=PhysicalUnit("electron_volt"),
        )
