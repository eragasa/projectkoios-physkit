"""Geometry and parameter binding for ordinary-cell ``sp3s*`` bonds."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.units import (
    PhysicalUnit,
    PintUnitConverter,
    UnitSystem,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterBondDirection:
    """Own one three-dimensional normalized Cartesian bond direction."""

    components: VectorQuantity

    def __post_init__(self) -> None:
        """Require one three-component vector of Euclidean norm one."""
        if type(self.components) is not VectorQuantity:
            raise TypeError("components must be VectorQuantity")
        if self.components.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("components must be nondimensional")
        if self.components.magnitude.shape != (3,):
            raise ValueError("components must have shape (3,)")
        tolerance = 128.0 * np.finfo(np.float64).eps
        if not math.isclose(
            float(np.linalg.norm(self.components.magnitude)),
            1.0,
            rel_tol=0.0,
            abs_tol=tolerance,
        ):
            raise ValueError("components must have unit Euclidean norm")


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterDirectedBond:
    """Identify one directed bond without embedding its Cartesian geometry."""

    ordinary_cell_displacement: LatticeDisplacement
    source_center_identifier: str
    target_center_identifier: str

    def __post_init__(self) -> None:
        """Validate one directed three-dimensional ordinary-cell bond."""
        if type(self.ordinary_cell_displacement) is not LatticeDisplacement:
            raise TypeError("ordinary_cell_displacement must be LatticeDisplacement")
        if self.ordinary_cell_displacement.dimension is not LatticeDimension.THREE:
            raise ValueError("sp3s* bonds require a three-dimensional lattice")
        for value, name in (
            (self.source_center_identifier, "source_center_identifier"),
            (self.target_center_identifier, "target_center_identifier"),
        ):
            if type(value) is not str:
                raise TypeError(f"{name} must be a built-in string")
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if (
            self.ordinary_cell_displacement.is_zero
            and self.source_center_identifier == self.target_center_identifier
        ):
            raise ValueError(
                "zero-cell self blocks belong to explicit onsite ownership"
            )

    def reversed(self) -> SlaterKosterDirectedBond:
        """Return the explicit Hermitian-reverse directed bond."""
        components = self.ordinary_cell_displacement.components
        if len(components) != 3:
            raise RuntimeError("validated three-dimensional displacement changed")
        return SlaterKosterDirectedBond(
            ordinary_cell_displacement=LatticeDisplacement(
                dimension=LatticeDimension.THREE,
                components=(-components[0], -components[1], -components[2]),
            ),
            source_center_identifier=self.target_center_identifier,
            target_center_identifier=self.source_center_identifier,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarBondGeometry:
    """Bind one directed bond to its physical Cartesian displacement.

    The directed bond's ordinary-cell displacement controls the Bloch phase and excludes
    orbital-center offsets. ``cartesian_displacement`` is the nonzero physical
    source-to-target vector and is the sole source of the normalized bond direction.
    """

    directed_bond: SlaterKosterDirectedBond
    cartesian_displacement: VectorQuantity
    unit_system: UnitSystem = UnitSystem.ATOMIC

    def __post_init__(self) -> None:
        """Validate the directed physical bond geometry."""
        if type(self.directed_bond) is not SlaterKosterDirectedBond:
            raise TypeError("directed_bond must be SlaterKosterDirectedBond")
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            raise ValueError("unit_system must provide physical length units")
        if type(self.cartesian_displacement) is not VectorQuantity:
            raise TypeError("cartesian_displacement must be VectorQuantity")
        if type(self.cartesian_displacement.unit) is not PhysicalUnit:
            raise ValueError("cartesian_displacement must carry a physical length unit")
        if not PintUnitConverter().compatible(
            left=self.cartesian_displacement.unit,
            right=PhysicalUnit(expression="meter"),
        ):
            raise ValueError("cartesian_displacement unit must represent length")
        if self.cartesian_displacement.unit != self.unit_system.length_unit:
            raise ValueError(
                "cartesian_displacement must use the selected unit-system length unit"
            )
        if self.cartesian_displacement.magnitude.shape != (3,):
            raise ValueError("cartesian_displacement must have shape (3,)")
        if float(np.linalg.norm(self.cartesian_displacement.magnitude)) == 0.0:
            raise ValueError("cartesian_displacement must be nonzero")

    @property
    def bond_direction(self) -> SlaterKosterBondDirection:
        """Return the immutable direction derived from the physical displacement."""
        displacement = self.cartesian_displacement.magnitude
        direction = displacement / np.linalg.norm(displacement)
        return SlaterKosterBondDirection(
            components=VectorQuantity(
                magnitude=np.asarray(direction, dtype=np.float64),
                unit=UnitSystem.NONDIMENSIONAL.length_unit,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarParameterizedBond:
    """Bind one bond geometry to an explicit shared parameter-group identifier."""

    geometry: HomopolarSp3sStarBondGeometry
    parameter_group_identifier: str

    def __post_init__(self) -> None:
        """Validate the concrete bond and parameter-group identifier."""
        if type(self.geometry) is not HomopolarSp3sStarBondGeometry:
            raise TypeError("geometry must be HomopolarSp3sStarBondGeometry")
        if type(self.parameter_group_identifier) is not str:
            raise TypeError("parameter_group_identifier must be a built-in string")
        if not self.parameter_group_identifier:
            raise ValueError("parameter_group_identifier must be nonempty")
