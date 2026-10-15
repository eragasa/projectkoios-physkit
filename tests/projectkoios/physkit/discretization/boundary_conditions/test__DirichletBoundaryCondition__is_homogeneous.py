"""Verification of homogeneous Dirichlet classification."""

from projectkoios.physkit.discretization.boundary_conditions import (
    DirichletBoundaryCondition,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity, Unitless


def test_exact_zero_is_homogeneous_and_retains_unit() -> None:
    boundary = DirichletBoundaryCondition(
        ScalarQuantity(0.0, PhysicalUnit("meter ** -0.5"))
    )

    assert boundary.is_homogeneous
    assert boundary.value.unit == PhysicalUnit("meter ** -0.5")
    assert not DirichletBoundaryCondition(
        ScalarQuantity(1.0, Unitless())
    ).is_homogeneous
