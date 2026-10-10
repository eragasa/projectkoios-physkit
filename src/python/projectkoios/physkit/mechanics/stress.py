"""Immutable unit-aware three-dimensional stress tensors."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from projectkoios.physkit.units import (
    MODEL_SYSTEM_UNIT_CONVERTER,
    MatrixQuantity,
    PhysicalUnit,
    ScalarQuantity,
)

_SI_STRESS_UNIT = PhysicalUnit("pascal")


class StressSignConvention(StrEnum):
    """Declare whether positive normal components denote tension or compression."""

    TENSION_POSITIVE = "tension-positive"
    COMPRESSION_POSITIVE = "compression-positive"


@dataclass(frozen=True, slots=True, eq=False)
class StressTensor:
    """Represent one immutable unit-aware Cartesian 3D stress tensor.

    Components retain their declared pressure-compatible unit and sign convention.
    Construction does not symmetrize or otherwise alter the supplied tensor.
    """

    components: MatrixQuantity
    sign_convention: StressSignConvention

    def __post_init__(self) -> None:
        if type(self.components) is not MatrixQuantity:
            raise TypeError("components must be a MatrixQuantity")
        if self.components.magnitude.shape != (3, 3):
            raise ValueError("stress tensor components must have shape (3, 3)")
        if type(self.components.unit) is not PhysicalUnit:
            raise TypeError("stress tensor components require a physical unit")
        if not MODEL_SYSTEM_UNIT_CONVERTER.compatible(
            self.components.unit,
            _SI_STRESS_UNIT,
        ):
            raise ValueError("stress tensor unit must have pressure dimensions")
        if type(self.sign_convention) is not StressSignConvention:
            raise TypeError("sign_convention must be a StressSignConvention")

    def to_unit(self, unit: PhysicalUnit) -> StressTensor:
        """Return an equivalent tensor expressed in a compatible pressure unit."""
        if type(unit) is not PhysicalUnit:
            raise TypeError("unit must be a PhysicalUnit")
        return StressTensor(
            components=MODEL_SYSTEM_UNIT_CONVERTER.convert_matrix(
                self.components,
                unit,
            ),
            sign_convention=self.sign_convention,
        )

    def to_si(self) -> StressTensor:
        """Return an equivalent tensor expressed in pascals."""
        return self.to_unit(_SI_STRESS_UNIT)

    def to_sign_convention(
        self,
        sign_convention: StressSignConvention,
    ) -> StressTensor:
        """Return an equivalent tensor using the requested explicit sign convention."""
        if type(sign_convention) is not StressSignConvention:
            raise TypeError("sign_convention must be a StressSignConvention")
        if sign_convention is self.sign_convention:
            return self
        return StressTensor(
            components=MatrixQuantity(
                -self.components.magnitude,
                self.components.unit,
            ),
            sign_convention=sign_convention,
        )

    @property
    def hydrostatic_pressure(self) -> ScalarQuantity:
        """Return mean normal pressure, positive under compression."""
        mean_normal_stress = float(np.trace(self.components.magnitude) / 3.0)
        if self.sign_convention is StressSignConvention.TENSION_POSITIVE:
            mean_normal_stress = -mean_normal_stress
        return ScalarQuantity(mean_normal_stress, self.components.unit)
