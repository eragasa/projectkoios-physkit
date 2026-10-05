"""Numerical verification for two-dimensional checkerboard Bloch bands."""

import numpy as np

from projectkoios.physkit.solidstate.bloch2d.checkerboard.model import (
    CheckerboardPotentialBloch2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.cell import (
    PeriodicPrimitiveCell2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.model import (
    PeriodicFreeParticle2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.tensor import (
    CartesianEffectiveMassTensor2D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    UnitSystem,
    VectorQuantity,
)


def _cell() -> PeriodicPrimitiveCell2D:
    return PeriodicPrimitiveCell2D(
        primitive_basis=MatrixQuantity(
            magnitude=10.0 * np.eye(2),
            unit=UnitSystem.METAL.length_unit,
        ),
        unit_system=UnitSystem.METAL,
    )


def _model(*, amplitude: float) -> CheckerboardPotentialBloch2D:
    return CheckerboardPotentialBloch2D(
        cell=_cell(),
        amplitude=ScalarQuantity(
            amplitude,
            UnitSystem.METAL.energy_unit,
        ),
    )


def _wave_vectors(values: np.ndarray) -> MatrixQuantity:
    return MatrixQuantity(
        magnitude=values,
        unit=_cell().reciprocal_basis.unit,
    )


def test__solve_bands__zero_amplitude_matches_free_particle_modes() -> None:
    model = _model(amplitude=0.0)
    wave_vector = np.array([[0.11, -0.07]])
    mode_limit = 2
    result = model.solve_bands(
        bloch_wave_vectors=_wave_vectors(wave_vector),
        reciprocal_mode_limit=mode_limit,
        band_count=8,
    )
    mode_indices = tuple(
        (first, second)
        for first in range(-mode_limit, mode_limit + 1)
        for second in range(-mode_limit, mode_limit + 1)
    )
    electron_mass = model.electron_mass
    free_model = PeriodicFreeParticle2D(
        cell=_cell(),
        effective_mass=CartesianEffectiveMassTensor2D(
            tensor=MatrixQuantity(
                magnitude=electron_mass.magnitude * np.eye(2),
                unit=electron_mass.unit,
            ),
            unit_system=UnitSystem.METAL,
        ),
    )
    free_energies = free_model.evaluate_modes(
        mode_indices=mode_indices,
        bloch_wave_vector=VectorQuantity(
            magnitude=wave_vector[0],
            unit=_cell().reciprocal_basis.unit,
        ),
    )

    np.testing.assert_allclose(
        result.energies.magnitude[0],
        np.sort(free_energies.energies.magnitude)[:8],
    )


def test__solve_bands__amplitude_sign_preserves_spectrum() -> None:
    wave_vectors = _wave_vectors(
        np.array(
            [
                [-0.2, 0.13],
                [0.0, 0.0],
                [0.17, -0.19],
            ]
        )
    )
    positive = _model(amplitude=0.7).solve_bands(
        bloch_wave_vectors=wave_vectors,
        reciprocal_mode_limit=2,
        band_count=9,
    )
    negative = _model(amplitude=-0.7).solve_bands(
        bloch_wave_vectors=wave_vectors,
        reciprocal_mode_limit=2,
        band_count=9,
    )

    np.testing.assert_allclose(
        positive.energies.magnitude,
        negative.energies.magnitude,
        rtol=2e-13,
        atol=2e-13,
    )
