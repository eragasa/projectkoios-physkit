"""Shared validation and numerics for separable cosine potentials."""

from __future__ import annotations

from itertools import product

import numpy as np

from projectkoios.physkit.numerics.typing.numpy.arrays import (
    ComplexVector,
    IntegerMatrix,
    RealMatrix,
    RealVector,
)
from projectkoios.physkit.units import (
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
)

_ENERGY_UNIT = PhysicalUnit("joule")


def validate_amplitude(amplitude: ScalarQuantity) -> None:
    """Require a physical energy or an explicitly nondimensional amplitude."""
    if not isinstance(amplitude, ScalarQuantity):
        raise TypeError("amplitude must be ScalarQuantity")
    if isinstance(amplitude.unit, PhysicalUnit) and not (
        ModelSystemUnitConverter.compatible(amplitude.unit, _ENERGY_UNIT)
    ):
        raise ValueError("amplitude must use units compatible with energy")
    if not isinstance(amplitude.unit, PhysicalUnit | Unitless):
        raise TypeError("amplitude must use PhysicalUnit or Unitless")


def validate_harmonics(harmonics: tuple[int, ...], dimension: int) -> None:
    """Require one positive built-in integer harmonic per dimension."""
    if type(harmonics) is not tuple:
        raise TypeError("harmonics must be a tuple")
    if len(harmonics) != dimension:
        raise ValueError(f"harmonics must contain {dimension} entries")
    if any(type(harmonic) is not int for harmonic in harmonics):
        raise TypeError("harmonics must contain built-in integers")
    if any(harmonic < 1 for harmonic in harmonics):
        raise ValueError("harmonics must be positive")


def validate_phases(phases: tuple[float, ...], dimension: int) -> None:
    """Require one finite built-in floating-point phase per dimension."""
    if type(phases) is not tuple:
        raise TypeError("phases must be a tuple")
    if len(phases) != dimension:
        raise ValueError(f"phases must contain {dimension} entries")
    if any(type(phase) is not float for phase in phases):
        raise TypeError("phases must contain built-in floats")
    if not np.all(np.isfinite(np.asarray(phases, dtype=np.float64))):
        raise ValueError("phases must be finite")


def cosine_product_fourier_support(
    harmonics: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    """Return the exact reciprocal-mode support of a cosine product."""
    return tuple(
        tuple(sign * harmonic for sign, harmonic in zip(signs, harmonics, strict=True))
        for signs in product((-1, 1), repeat=len(harmonics))
    )


def evaluate_cosine_product(
    *,
    amplitude: float,
    harmonics: tuple[int, ...],
    phases: tuple[float, ...],
    fractional_coordinates: RealMatrix,
) -> RealVector:
    r"""Evaluate $A\prod_j\cos(2\pi h_j u_j+\phi_j)$."""
    arguments = (
        2.0
        * np.pi
        * fractional_coordinates
        * np.asarray(harmonics, dtype=np.float64)[np.newaxis, :]
        + np.asarray(phases, dtype=np.float64)[np.newaxis, :]
    )
    return np.asarray(
        amplitude * np.prod(np.cos(arguments), axis=1),
        dtype=np.float64,
    )


def cosine_product_fourier_coefficients(
    *,
    amplitude: float,
    harmonics: tuple[int, ...],
    phases: tuple[float, ...],
    mode_indices: IntegerMatrix,
) -> ComplexVector:
    r"""Return coefficients in $V(u)=\sum_m V_m e^{i2\pi m\cdot u}$."""
    harmonic_array = np.asarray(harmonics, dtype=np.int64)
    coupled = np.all(np.abs(mode_indices) == harmonic_array[np.newaxis, :], axis=1)
    coefficients = np.zeros(mode_indices.shape[0], dtype=np.complex128)
    signs = mode_indices[coupled] / harmonic_array[np.newaxis, :]
    phase_array = np.asarray(phases, dtype=np.float64)
    coefficients[coupled] = (
        amplitude / (2.0 ** len(harmonics)) * np.exp(1j * (signs @ phase_array))
    )
    return coefficients
