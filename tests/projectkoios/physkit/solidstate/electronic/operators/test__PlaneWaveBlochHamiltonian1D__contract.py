"""Contract tests for the one-dimensional plane-wave Hamiltonian."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis1D
from projectkoios.physkit.periodic.lattice import (
    DirectLattice1D,
    ReciprocalLattice1D,
)
from projectkoios.physkit.solidstate.electronic.models import (
    CosinePotentialFourier1D,
    PlaneWaveBlochModel1D,
)
from projectkoios.physkit.solidstate.electronic.operators import (
    PlaneWaveBlochHamiltonian1D,
)


def test__construct_matrix__uses_canonical_reciprocal_vectors() -> None:
    """Kinetic and potential terms correlate to one canonical mode basis."""
    direct = DirectLattice1D(2.0 * np.pi)
    reciprocal = ReciprocalLattice1D.from_direct_lattice(direct)
    basis = ReciprocalModeBasis1D.from_mode_limit(reciprocal, 1)
    model = PlaneWaveBlochModel1D(
        lattice=direct,
        potential=CosinePotentialFourier1D(amplitude=2.0, lattice=direct),
    )

    matrix = PlaneWaveBlochHamiltonian1D(model, basis).construct_matrix(0.0)

    np.testing.assert_allclose(
        matrix,
        np.asarray(
            [
                [0.5, 1.0, 0.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 0.5],
            ],
            dtype=np.complex128,
        ),
        atol=1.0e-15,
    )
