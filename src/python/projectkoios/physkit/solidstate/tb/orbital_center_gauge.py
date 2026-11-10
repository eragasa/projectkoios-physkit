r"""Gauge transformation for explicit orbital-center embeddings.

For ordered fractional orbital centers :math:`\tau_a` and fractional reciprocal
coordinate :math:`q`, this module constructs

.. math::

   D_{aa}(q)=\exp(+2\pi i\,q\cdot\tau_a).

A cell-gauge Bloch matrix transforms to the full-center representation as
:math:`H_{\mathrm{full}}=D^\dagger H_{\mathrm{cell}}D`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.core.results import ResultsObject
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    MatrixQuantity,
    UnitSystem,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CellToFullCenterGaugeTransformationDefinition:
    """Bind ordered fractional orbital centers and one reciprocal coordinate."""

    orbital_centers_fractional: MatrixQuantity
    fractional_wave_vector: VectorQuantity

    def __post_init__(self) -> None:
        """Validate the retained fractional-coordinate quantities."""
        if type(self.orbital_centers_fractional) is not MatrixQuantity:
            raise TypeError("orbital_centers_fractional must be MatrixQuantity")
        if (
            self.orbital_centers_fractional.unit
            != UnitSystem.NONDIMENSIONAL.length_unit
        ):
            raise ValueError("orbital_centers_fractional must be nondimensional")
        if self.orbital_centers_fractional.magnitude.ndim != 2:
            raise ValueError("orbital_centers_fractional must be a matrix")
        if self.orbital_centers_fractional.magnitude.shape[1] != 3:
            raise ValueError(
                "orbital_centers_fractional must have shape (basis_size, 3)"
            )
        if type(self.fractional_wave_vector) is not VectorQuantity:
            raise TypeError("fractional_wave_vector must be VectorQuantity")
        if self.fractional_wave_vector.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("fractional_wave_vector must be nondimensional")
        if self.fractional_wave_vector.magnitude.shape != (3,):
            raise ValueError("fractional_wave_vector must have shape (3,)")


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class CellToFullCenterGaugeTransformation(ResultsObject):
    """Correlate one gauge definition with its diagonal transformation matrix."""

    definition: CellToFullCenterGaugeTransformationDefinition
    matrix: ComplexMatrixQuantity

    def __post_init__(self) -> None:
        """Validate the represented square transformation."""
        if type(self.definition) is not CellToFullCenterGaugeTransformationDefinition:
            raise TypeError(
                "definition must be CellToFullCenterGaugeTransformationDefinition"
            )
        if type(self.matrix) is not ComplexMatrixQuantity:
            raise TypeError("matrix must be ComplexMatrixQuantity")
        if self.matrix.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("matrix must be nondimensional")
        basis_size = self.definition.orbital_centers_fractional.magnitude.shape[0]
        if self.matrix.magnitude.shape != (basis_size, basis_size):
            raise ValueError("matrix shape must match the ordered orbital centers")


class CellToFullCenterGaugeTransformationConstructor:
    """Construct one explicit diagonal cell-to-full-center gauge transformation."""

    __slots__ = ()

    def execute(
        self,
        *,
        definition: CellToFullCenterGaugeTransformationDefinition,
    ) -> CellToFullCenterGaugeTransformation:
        r"""Return :math:`D_{aa}(q)=\exp(+2\pi i q\cdot\tau_a)`."""
        if type(definition) is not CellToFullCenterGaugeTransformationDefinition:
            raise TypeError(
                "definition must be CellToFullCenterGaugeTransformationDefinition"
            )
        phases = np.exp(
            2j
            * np.pi
            * (
                definition.orbital_centers_fractional.magnitude
                @ definition.fractional_wave_vector.magnitude
            )
        )
        return CellToFullCenterGaugeTransformation(
            definition=definition,
            matrix=ComplexMatrixQuantity(
                magnitude=np.diag(phases),
                unit=UnitSystem.NONDIMENSIONAL.length_unit,
            ),
        )
