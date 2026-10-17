"""Contract tests for the one-dimensional plane-wave band solver."""

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.basis import ReciprocalModeBasis1D
from projectkoios.physkit.periodic.kpoints import KPointPath1D
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
from projectkoios.physkit.solidstate.electronic.solvers import PlaneWaveBlochSolver1D


def test__solve__retains_grid_and_canonical_reciprocal_vectors() -> None:
    """The solver returns complete evidence correlated to its canonical basis."""
    direct = DirectLattice1D(2.0 * np.pi)
    basis = ReciprocalModeBasis1D.from_mode_limit(
        ReciprocalLattice1D.from_direct_lattice(direct),
        1,
    )
    model = PlaneWaveBlochModel1D(
        lattice=direct,
        potential=CosinePotentialFourier1D(amplitude=0.5, lattice=direct),
    )
    hamiltonian = PlaneWaveBlochHamiltonian1D(model, basis)
    solver = PlaneWaveBlochSolver1D(hamiltonian)
    path = KPointPath1D(np.asarray([-0.25, 0.0, 0.25]))

    result = solver.solve(path, number_of_bands=2)

    np.testing.assert_array_equal(result.k_points, path.values)
    np.testing.assert_array_equal(result.reciprocal_vectors, basis.vectors)
    assert result.energies.shape == (3, 2)
    assert result.eigenvectors.shape == (3, 3, 2)
    for point_index, wavevector in enumerate(path.values):
        expected = np.linalg.eigvalsh(hamiltonian.construct_matrix(wavevector))[:2]
        np.testing.assert_allclose(result.energies[point_index], expected)
