"""Fourier-coefficient tests for ``CosinePotential1D``."""

import numpy as np

from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.units import PhysicalUnit, ScalarQuantity


def test__fourier_coefficients__returns_only_conjugate_harmonic_pair() -> None:
    """A shifted real cosine has the expected conjugate coefficients."""
    unit = PhysicalUnit("electron_volt")
    potential = CosinePotential1D(
        amplitude=ScalarQuantity(4.0, unit),
        harmonic=2,
        phase=np.pi / 3.0,
    )

    evaluation = potential.fourier_coefficients(
        mode_indices=((-2,), (-1,), (0,), (1,), (2,))
    )

    expected = np.asarray(
        [2.0 * np.exp(-1j * np.pi / 3.0), 0.0, 0.0, 0.0,
         2.0 * np.exp(1j * np.pi / 3.0)],
        dtype=np.complex128,
    )
    assert evaluation.coefficients.unit == unit
    np.testing.assert_allclose(evaluation.coefficients.magnitude, expected)
