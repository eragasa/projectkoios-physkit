"""Represented checkerboard Fourier-coupling verification."""

import numpy as np

from projectkoios.physkit.solidstate.bloch2d.cell import BlochCell2D
from projectkoios.physkit.solidstate.bloch2d.checkerboard.model import (
    CheckerboardPotentialBloch2D,
)
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    ScalarQuantity,
    Unitless,
    UnitSystem,
    VectorQuantity,
)


def test__construct_hamiltonian__has_checkerboard_fourier_couplings() -> None:
    amplitude = 0.8
    cell = BlochCell2D(
        direct_basis=MatrixQuantity(
            magnitude=10.0 * np.eye(2),
            unit=UnitSystem.METAL.length_unit,
        ),
        unit_system=UnitSystem.METAL,
    )
    model = CheckerboardPotentialBloch2D(
        cell=cell,
        amplitude=ScalarQuantity(amplitude, UnitSystem.METAL.energy_unit),
    )
    represented = model.construct_hamiltonian(
        bloch_wave_vector=VectorQuantity(
            magnitude=np.zeros(2),
            unit=cell.reciprocal_basis.unit,
        ),
        reciprocal_mode_limit=1,
    )
    assert isinstance(
        represented.represented_hamiltonian,
        ComplexSparseMatrixQuantity,
    )
    assert (
        represented.represented_hamiltonian.nonzero_count
        < len(represented.mode_indices) ** 2
    )
    matrix = represented.represented_hamiltonian.to_csr().toarray()
    modes = represented.mode_indices
    zero_index = modes.index((0, 0))

    for index, mode in enumerate(modes):
        if abs(mode[0]) == 1 and abs(mode[1]) == 1:
            assert matrix[zero_index, index] == amplitude / 4.0
        elif index != zero_index:
            assert matrix[zero_index, index] == 0.0
    np.testing.assert_allclose(matrix, matrix.conjugate().T)


def test__construct_hamiltonian__uses_unitless_hbar_and_electron_mass_one() -> None:
    cell = BlochCell2D(
        direct_basis=MatrixQuantity(
            magnitude=np.eye(2),
            unit=Unitless(),
        ),
        unit_system=UnitSystem.NONDIMENSIONAL,
    )
    model = CheckerboardPotentialBloch2D(
        cell=cell,
        amplitude=ScalarQuantity(0.0, Unitless()),
    )

    represented = model.construct_hamiltonian(
        bloch_wave_vector=VectorQuantity(
            magnitude=np.zeros(2),
            unit=Unitless(),
        ),
        reciprocal_mode_limit=1,
    )

    mode_matrix = np.asarray(represented.mode_indices, dtype=np.float64)
    reciprocal_vectors = mode_matrix @ cell.reciprocal_basis.magnitude.T
    expected_energies = 0.5 * np.einsum(
        "ni,ni->n",
        reciprocal_vectors,
        reciprocal_vectors,
    )
    np.testing.assert_allclose(
        represented.represented_hamiltonian.to_csr().diagonal(),
        expected_energies,
    )
    assert represented.represented_hamiltonian.unit == Unitless()
    assert model.electron_mass == ScalarQuantity(1.0, Unitless())
