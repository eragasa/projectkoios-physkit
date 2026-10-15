# physkit/qm/solver1d.py
# Eugene Ragasa, 2026

from __future__ import annotations

import numpy as np
import scipy.sparse.linalg as spla
from numpy.typing import NDArray
from scipy.sparse import csr_matrix

from projectkoios.physkit.core.boundaries import (
    BoundaryCondition,
)
from projectkoios.physkit.core.boundaries import (
    DirichletBoundaryCondition as DirichletBC,
)
from projectkoios.physkit.core.grids import CartesianGrid1D
from projectkoios.physkit.numerics.finite_difference import Laplacian1D
from projectkoios.physkit.qm.potentials import Potential1D

FloatArray = NDArray[np.float64]
_DEFAULT_DIRICHLET_BC = DirichletBC(value=0.0)


class TISEAnalytical1D:
    def __init__(self, x_lower: float, x_upper: float, mass: float, hbar: float):
        self.x_min: float = x_lower
        self.x_max: float = x_upper
        self.mass: float = mass
        self.hbar: float = mass


class TISESymbolic1D:
    def __init__(self): ...


class TISESolver1D:
    """
    One-dimensional finite-difference Schrödinger solver.

    This class solves the time-independent Schrödinger equation

        H psi = E psi

    on a CartesianGrid1D.

    The Hamiltonian is

        H = T + V

    where

        T = -(hbar^2 / 2m) d^2/dx^2

    and V is represented as a diagonal matrix on the active grid.

    First implementation
    --------------------
    This solver currently assumes that Laplacian1D supports only homogeneous
    Dirichlet boundary conditions:

        psi(x_min) = 0
        psi(x_max) = 0

    Therefore the active state vector excludes the boundary grid points.
    """

    def __init__(
        self,
        grid: CartesianGrid1D,
        potential: Potential1D,
        *,
        mass: float = 1.0,
        hbar: float = 1.0,
        bc_x_min: BoundaryCondition = _DEFAULT_DIRICHLET_BC,
        bc_x_max: BoundaryCondition = _DEFAULT_DIRICHLET_BC,
    ) -> None:
        self.mass: float = mass
        self.hbar: float = hbar
        self.grid: CartesianGrid1D = grid
        self.potential: Potential1D = potential
        self.bc_x_min: BoundaryCondition = bc_x_min
        self.bc_x_max: BoundaryCondition = bc_x_max

        self.check_args()

        self.laplacian: Laplacian1D = Laplacian1D(
            grid=self.grid,
            bc_x_min=self.bc_x_min,
            bc_x_max=self.bc_x_max,
        )

        self.laplacian_matrix: csr_matrix = self.laplacian.build()
        self.kinetic_matrix: csr_matrix = self.build_kinetic_matrix()
        self.potential_matrix: csr_matrix = self.build_potential_matrix()
        self.hamiltonian_matrix: csr_matrix = self.build_hamiltonian_matrix()

        self.eigval: np.ndarray | None = None
        self.eigvec: np.ndarray | None = None
        self.eig_n_states: int | None = None

        self.energies: np.ndarray | None = None
        self.psi: np.ndarray | None = None

    @property
    def x_active(self) -> np.ndarray:
        return self.laplacian.x_active

    def check_potential(self) -> None:
        if not isinstance(self.potential, Potential1D):
            raise TypeError("potential must be a Potential1D")

        # check the grid
        V_grid = self.potential.eval(self.grid.x)
        if not isinstance(V_grid, np.ndarray):
            raise TypeError("potential.eval(grid.x) must return a NumPy array.")
        if V_grid.shape != self.grid.shape:
            raise ValueError("potential.eval(grid.x) must have shape grid.shape.")
        if not np.all(np.isfinite(V_grid)):
            raise ValueError("potential values must be finite.")

    def check_mass(self) -> None:
        if not np.isfinite(self.mass):
            raise ValueError("mass must be finite.")

        if self.mass <= 0.0:
            raise ValueError("mass must be positive.")

        if not np.isfinite(self.hbar):
            raise ValueError("hbar must be finite.")

    def check_hbar(self) -> None:
        if not np.isfinite(self.hbar):
            raise ValueError("hbar must be finite.")
        if self.hbar <= 0.0:
            raise ValueError("hbar must be positive.")

    def check_args(self) -> None:
        if not isinstance(self.grid, CartesianGrid1D):
            raise TypeError("grid must be a CartesianGrid1D.")

        self.check_potential()
        self.check_mass()
        self.check_hbar()

        if not isinstance(self.bc_x_min, BoundaryCondition):
            raise TypeError("bc_x_min must be a BoundaryCondition.")

        if not isinstance(self.bc_x_max, BoundaryCondition):
            raise TypeError("bc_x_max must be a BoundaryCondition.")

    def build_kinetic_matrix(self) -> csr_matrix:
        """
        Build the kinetic-energy matrix.

        The continuous kinetic-energy operator is

            T = -(hbar^2 / 2m) d^2/dx^2
        """
        return -(self.hbar**2 / (2.0 * self.mass)) * self.L

    def build_potential_matrix(self) -> csr_matrix:
        """
        Build the diagonal potential-energy matrix.
        """
        return self.potential.sparse_matrix(self.x_active)

    def build_hamiltonian_matrix(self) -> csr_matrix:
        """
        Build the Hamiltonian matrix.

            H = T + V
        """
        return self.T + self.V

    def solve(
        self,
        *,
        eig_n_states: int = 6,
    ) -> tuple[FloatArray, FloatArray]:
        """
        Solve the stationary Schrödinger eigenvalue problem.

        Parameters
        ----------
        n_states:
            Number of lowest-energy states to compute.

        Returns
        -------
        energies:
            Eigenvalues sorted from lowest to highest.
        states:
            Eigenvectors stored columnwise. states[:, n] is eigenstate n.
        """
        self.check_eig_n_states(eig_n_states=eig_n_states)
        self.eig_n_states = eig_n_states

        self.diagonalize()

        idx = np.argsort(self.eigval)
        self.energies = self.eigval[idx]
        self.psi = self.eigvec[:, idx]

        self.normalize_states()

        return self.energies, self.psi

    def diagonalize(self) -> tuple[np.ndarray, np.ndarray]:
        self.eigval, self.eigvec = spla.eigsh(
            self.H,
            k=self.eig_n_states,
            which="SA",
        )
        return self.eigval, self.eigvec

    def check_eig_n_states(self, eig_n_states: int):
        if not isinstance(eig_n_states, int):
            raise TypeError("eig_n_states must be an integer.")

        if eig_n_states < 1:
            raise ValueError("eig_n_states must be at least 1.")

        if eig_n_states >= self.H.shape[0]:
            raise ValueError("eig_n_states must be smaller than the matrix size.")

    def normalize_states(self) -> FloatArray:
        """
        Normalize eigenstates using the discrete integral convention.

        The normalization condition is

            integral |psi(x)|^2 dx = 1

        represented numerically by

            sum_i |psi_i|^2 dx = 1
        """
        if not isinstance(self.psi, np.ndarray):
            raise TypeError("states must be a NumPy array.")

        norms = []
        psi_norm = self.eigvec.copy()

        for n in range(psi_norm.shape[1]):
            norm_sq = np.vdot(psi_norm[:, n], psi_norm[:, n]) * self.grid.dx

            if not np.isclose(np.imag(norm_sq), 0.0):
                raise ValueError("State norm has a nonzero imaginary component.")

            norm = np.sqrt(np.real(norm_sq))
            norms.append(norm)

            if norm == 0.0:
                raise ValueError("Cannot normalize a zero eigenstate.")

            psi_norm[:, n] = psi_norm[:, n] / norm

        self.norms = norms
        self.psi = psi_norm

        return norms, psi_norm
