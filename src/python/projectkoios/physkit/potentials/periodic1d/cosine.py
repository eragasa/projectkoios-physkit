r"""One-dimensional periodic cosine potential."""

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
class CosinePotential1D(BasePeriodicPotential):
    r"""Represent $V(u)=A\cos(2\pi h u+\phi)$ on a unit cell."""

    amplitude: ScalarQuantity
    harmonic: int = 1
    phase: float = 0.0

    def __post_init__(self) -> None:
        """Check amplitude, harmonic, and phase."""
        validate_amplitude(self.amplitude)
        validate_harmonics((self.harmonic,), self.dimension)
        validate_phases((self.phase,), self.dimension)

    @property
    def dimension(self) -> int:
        """Return the one-dimensional coordinate count."""
        return 1

    @property
    def energy_unit(self) -> ModelSystemUnit:
        """Return the amplitude unit."""
        return self.amplitude.unit

    @property
    def fourier_support(self) -> tuple[tuple[int, ...], ...]:
        """Return the two coupled reciprocal harmonics."""
        return cosine_product_fourier_support((self.harmonic,))

    def _evaluate_magnitudes(
        self,
        fractional_coordinates: RealMatrix,
    ) -> RealVector:
        return evaluate_cosine_product(
            amplitude=self.amplitude.magnitude,
            harmonics=(self.harmonic,),
            phases=(self.phase,),
            fractional_coordinates=fractional_coordinates,
        )

    def _fourier_coefficient_magnitudes(
        self,
        mode_indices: IntegerMatrix,
    ) -> ComplexVector:
        return cosine_product_fourier_coefficients(
            amplitude=self.amplitude.magnitude,
            harmonics=(self.harmonic,),
            phases=(self.phase,),
            mode_indices=mode_indices,
        )
