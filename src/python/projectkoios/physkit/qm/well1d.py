# physkit/qm/well1d.py

from dataclasses import dataclass

import numpy as np

from ..core.operator import LinearOperator1D
from ..core.state import Grid1D, Wavefunction1D


@dataclass
class InfiniteSquareWell1D:
    """
    Physics container for the 1D infinite square well on (0, L).
    Builds a Grid1D and the corresponding Hamiltonian operator.
    """

    L: float
    n_points: int
    m: float = 1.0
    hbar: float = 1.0

    def make_grid(self) -> Grid1D:
        """
        Create a Grid1D for internal points only (Dirichlet at 0 and L).
        """
        # total points including boundaries: N_total
        N_total = self.n_points
        x_full = np.linspace(0.0, self.L, N_total)
        # internal points: exclude boundaries
        x_internal = x_full[1:-1]
        return Grid1D(x_internal)

    def make_hamiltonian(self) -> InfiniteSquareWellHamiltonian1D:
        grid = self.make_grid()
        return InfiniteSquareWellHamiltonian1D(grid=grid, m=self.m, hbar=self.hbar)


class InfiniteSquareWellHamiltonian1D(LinearOperator1D):
    """
    Finite-difference Hamiltonian for the 1D infinite square well.

    Acts on Wavefunction1D defined on internal grid points (0, L) with
    implicit Dirichlet boundary conditions psi(0)=psi(L)=0.
    """

    def __init__(self, grid: Grid1D, m: float = 1.0, hbar: float = 1.0):
        super().__init__(grid)
        self.m = m
        self.hbar = hbar

    def _build_matrix(self) -> np.ndarray:
        x = self.grid.x
        N = x.size
        dx = self.grid.dx

        # Second derivative with Dirichlet BCs on full grid translates to
        # this tridiagonal stencil on internal points:
        # (psi_{i+1} - 2 psi_i + psi_{i-1}) / dx^2
        diag = np.full(N, -2.0)
        off = np.ones(N - 1)
        D2 = (np.diag(diag) + np.diag(off, k=1) + np.diag(off, k=-1)) / dx**2

        # Hamiltonian H = -(ħ^2 / 2m) D2
        factor = -(self.hbar**2) / (2.0 * self.m)
        H = factor * D2
        return H


def analytic_energy_levels(n, L: float, m: float = 1.0, hbar: float = 1.0):
    """
    Analytic levels for comparison:
    E_n = (ħ^2 π^2 n^2) / (2 m L^2)
    """
    n = np.asarray(n, dtype=float)
    return (hbar**2 * np.pi**2 * n**2) / (2.0 * m * L**2)


def reconstruct_with_boundaries(wf: Wavefunction1D, L: float):
    """
    For plotting: embed internal wavefunction into a full vector with
    psi(0)=psi(L)=0, and construct the corresponding full x-grid.
    """
    x_internal = wf.grid.x
    # infer N_total = N_internal + 2
    N_internal = x_internal.size
    N_total = N_internal + 2
    dx = L / (N_total - 1)  # noqa: F841 - Retain the derived grid spacing.
    x_full = np.linspace(0.0, L, N_total)

    psi_full = np.zeros(N_total, dtype=complex)
    psi_full[1:-1] = wf.values

    return x_full, psi_full
