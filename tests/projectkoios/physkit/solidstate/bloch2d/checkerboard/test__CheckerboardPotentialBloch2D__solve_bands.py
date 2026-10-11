"""Numerical verification for two-dimensional checkerboard Bloch bands."""

import numpy as np
import pytest

from projectkoios.physkit.solidstate.bloch2d.cell import BlochCell2D
from projectkoios.physkit.solidstate.bloch2d.checkerboard import model as model_module
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


def _cell() -> BlochCell2D:
    return BlochCell2D(
        direct_basis=MatrixQuantity(
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
        cell=PeriodicPrimitiveCell2D(
            primitive_basis=MatrixQuantity(
                magnitude=10.0 * np.eye(2),
                unit=UnitSystem.METAL.length_unit,
            ),
            unit_system=UnitSystem.METAL,
        ),
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


def test__solve_bands__uses_sparse_solver_for_a_partial_spectrum(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[int, int]] = []
    original_eigsh = model_module.sparse_linalg.eigsh

    def recording_eigsh(operator: object, *, k: int, **kwargs: object) -> np.ndarray:
        calls.append((operator.shape[0], k))  # type: ignore[attr-defined]
        return original_eigsh(operator, k=k, **kwargs)

    monkeypatch.setattr(model_module.sparse_linalg, "eigsh", recording_eigsh)

    result = _model(amplitude=0.4).solve_bands(
        bloch_wave_vectors=_wave_vectors(np.array([[0.03, -0.07]])),
        reciprocal_mode_limit=2,
        band_count=6,
    )

    assert result.energies.magnitude.shape == (1, 6)
    assert calls == [(25, 6)]


def test__solve_bands__sparse_partial_spectrum_matches_dense_reference() -> None:
    model = _model(amplitude=0.7)
    wave_vector = np.array([0.08, -0.04])
    represented = model.construct_hamiltonian(
        bloch_wave_vector=VectorQuantity(
            magnitude=wave_vector,
            unit=_cell().reciprocal_basis.unit,
        ),
        reciprocal_mode_limit=2,
    )
    dense_eigenvalues = np.linalg.eigvalsh(
        represented.represented_hamiltonian.to_csr().toarray()
    )

    result = model.solve_bands(
        bloch_wave_vectors=_wave_vectors(wave_vector[np.newaxis, :]),
        reciprocal_mode_limit=2,
        band_count=7,
    )

    np.testing.assert_allclose(
        result.energies.magnitude[0],
        dense_eigenvalues[:7],
        rtol=2e-13,
        atol=2e-13,
    )


def test__solve_bands__complete_spectrum_uses_explicit_dense_boundary() -> None:
    model = _model(amplitude=0.5)
    wave_vector = np.array([0.05, 0.02])
    represented = model.construct_hamiltonian(
        bloch_wave_vector=VectorQuantity(
            magnitude=wave_vector,
            unit=_cell().reciprocal_basis.unit,
        ),
        reciprocal_mode_limit=1,
    )
    expected = np.linalg.eigvalsh(
        represented.represented_hamiltonian.to_csr().toarray()
    )

    result = model.solve_bands(
        bloch_wave_vectors=_wave_vectors(wave_vector[np.newaxis, :]),
        reciprocal_mode_limit=1,
        band_count=9,
    )

    np.testing.assert_allclose(result.energies.magnitude[0], expected)


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
