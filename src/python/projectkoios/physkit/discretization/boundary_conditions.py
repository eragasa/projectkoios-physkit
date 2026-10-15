"""Unit-aware boundary conditions migrated from ksdft2effmass model systems."""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from projectkoios.physkit.core.boundaries import (
    BoundaryCondition,
    BoundaryConditionType,
)
from projectkoios.physkit.units import ScalarQuantity


@dataclass(frozen=True, slots=True)
class DirichletBoundaryCondition(BoundaryCondition):
    """Represent one constant value prescribed at selected boundary points."""

    value: ScalarQuantity
    type: ClassVar[BoundaryConditionType] = BoundaryConditionType.DIRICHLET

    def __post_init__(self) -> None:
        if not isinstance(self.value, ScalarQuantity):
            raise TypeError("value must be ScalarQuantity")

    @property
    def condition_kind(self) -> str:
        """Return the exact boundary-condition kind identifier."""
        return "dirichlet"

    @property
    def is_homogeneous(self) -> bool:
        """Return whether the prescribed boundary value is exactly zero."""
        return self.value.magnitude == 0.0
