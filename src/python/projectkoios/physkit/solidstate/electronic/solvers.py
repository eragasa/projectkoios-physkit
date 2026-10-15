# src/python/projectkoios/physkit/solidstate/electronic/solvers.py

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.kpoints import KPointPath1D
from projectkoios.physkit.solidstate.electronic.operators import (
    PlaneWaveBlochHamiltonian1D,
)
from projectkoios.physkit.solidstate.electronic.results import (
    ElectronicBandStructureResult,
)

type FloatArray = NDArray[np.float64]
type ComplexArray = NDArray[np.complex128]


class PlaneWaveBlochSolver1D:
    """
    Eigensolver for a one-dimensional electronic Bloch Hamiltonian.

    Parameters
    ----------
    hamiltonian:
        Plane-wave Bloch Hamiltonian builder.
    """

    def __init__(
        self,
        hamiltonian: PlaneWaveBlochHamiltonian1D,
    ) -> None:
        self.hamiltonian = hamiltonian

    def solve_at_wavevector(
        self,
        bloch_wavevector: float,
        *,
        number_of_bands: int,
    ) -> tuple[FloatArray, ComplexArray]:
        """
        Solve the electronic eigenproblem at one wavevector.

        Parameters
        ----------
        bloch_wavevector:
            Bloch wavevector ``k``.
        number_of_bands:
            Number of lowest-energy bands retained.

        Returns
        -------
        energies:
            Lowest band energies.
        eigenvectors:
            Corresponding plane-wave coefficient vectors.
        """
        if not 1 <= number_of_bands <= self.hamiltonian.basis.size:
            raise ValueError("number_of_bands must be between 1 and the basis size.")

        matrix = self.hamiltonian.construct_matrix(bloch_wavevector)

        energies, eigenvectors = np.linalg.eigh(matrix)

        return (
            energies[:number_of_bands],
            eigenvectors[:, :number_of_bands],
        )

    def solve(
        self,
        k_path: KPointPath1D,
        *,
        number_of_bands: int,
    ) -> ElectronicBandStructureResult:
        """
        Solve the electronic band structure over a wavevector path.

        Parameters
        ----------
        k_path:
            Ordered Bloch-wavevector path.
        number_of_bands:
            Number of bands retained at every wavevector.

        Returns
        -------
        ElectronicBandStructureResult
            Complete electronic band-structure result.
        """
        number_of_k_points = k_path.size
        basis_size = self.hamiltonian.basis.size

        energies = np.zeros(
            (
                number_of_k_points,
                number_of_bands,
            ),
            dtype=np.float64,
        )

        eigenvectors = np.zeros(
            (
                number_of_k_points,
                basis_size,
                number_of_bands,
            ),
            dtype=np.complex128,
        )

        for k_index, bloch_wavevector in enumerate(k_path.values):
            (
                energies[k_index],
                eigenvectors[k_index],
            ) = self.solve_at_wavevector(
                bloch_wavevector,
                number_of_bands=number_of_bands,
            )

        return ElectronicBandStructureResult(
            k_points=k_path.values.copy(),
            energies=energies,
            eigenvectors=eigenvectors,
            reciprocal_vectors=(self.hamiltonian.basis.vectors.copy()),
        )
