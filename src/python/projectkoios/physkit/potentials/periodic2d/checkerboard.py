r"""Two-dimensional separable checkerboard potential."""

from __future__ import annotations

from dataclasses import dataclass

from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexVector,
    IntegerMatrix,
    RealMatrix,
    RealVector,
)
from projectkoios.physkit.potentials.periodic._cosine_product import (
    cosine_product_fourier_coefficients,
    cosine_product_fourier_support,
    evaluate_cosine_product,
    validate_amplitude,
    validate_harmonics,
    validate_phases,
)
from projectkoios.physkit.potentials.periodic.base import BasePeriodicPotential
from projectkoios.physkit.units import ModelSystemUnit, ScalarQuantity


@dataclass(frozen=True, slots=True, kw_only=True)
class CheckerboardPotential2D(BasePeriodicPotential):
    r"""Represent $V(u)=A\prod_{j=1}^2\cos(2\pi h_j u_j+\phi_j)$."""

    amplitude: ScalarQuantity
    harmonics: tuple[int, int] = (1, 1)
    phases: tuple[float, float] = (0.0, 0.0)

    def __post_init__(self) -> None:
        """Check amplitude, harmonics, and phases."""
        validate_amplitude(self.amplitude)
        validate_harmonics(self.harmonics, self.dimension)
        validate_phases(self.phases, self.dimension)

    @property
    def dimension(self) -> int:
        """Return the two-dimensional coordinate count."""
        return 2

    @property
    def energy_unit(self) -> ModelSystemUnit:
        """Return the amplitude unit."""
        return self.amplitude.unit

    @property
    def fourier_support(self) -> tuple[tuple[int, ...], ...]:
        """Return the four coupled reciprocal modes."""
        return cosine_product_fourier_support(self.harmonics)

    def _evaluate_magnitudes(
        self,
        fractional_coordinates: RealMatrix,
    ) -> RealVector:
        return evaluate_cosine_product(
            amplitude=self.amplitude.magnitude,
            harmonics=self.harmonics,
            phases=self.phases,
            fractional_coordinates=fractional_coordinates,
        )

    def _fourier_coefficient_magnitudes(
        self,
        mode_indices: IntegerMatrix,
    ) -> ComplexVector:
        return cosine_product_fourier_coefficients(
            amplitude=self.amplitude.magnitude,
            harmonics=self.harmonics,
            phases=self.phases,
            mode_indices=mode_indices,
        )
