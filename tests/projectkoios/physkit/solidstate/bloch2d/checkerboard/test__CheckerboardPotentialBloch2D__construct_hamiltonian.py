"""Represented checkerboard Fourier-coupling verification."""

import numpy as np

from projectkoios.physkit.solidstate.bloch2d.checkerboard.model import (
    CheckerboardPotentialBloch2D,
)
from projectkoios.physkit.solidstate.semiconductors.effmass_2d.cell import (
    PeriodicPrimitiveCell2D,
)
from projectkoios.physkit.units import (
    MatrixQuantity,
    ScalarQuantity,
    UnitSystem,
    VectorQuantity,
)


def test__construct_hamiltonian__has_checkerboard_fourier_couplings() -> None:
    amplitude = 0.8
    cell = PeriodicPrimitiveCell2D(
        primitive_basis=MatrixQuantity(
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
    matrix = represented.represented_hamiltonian.magnitude
    modes = represented.mode_indices
    zero_index = modes.index((0, 0))

    for index, mode in enumerate(modes):
        if abs(mode[0]) == 1 and abs(mode[1]) == 1:
            assert matrix[zero_index, index] == amplitude / 4.0
        elif index != zero_index:
            assert matrix[zero_index, index] == 0.0
    np.testing.assert_allclose(matrix, matrix.conjugate().T)
