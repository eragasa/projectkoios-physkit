# physkit/core/state.py

from dataclasses import dataclass

import numpy as np

# physkit/core/state.py


@dataclass
class Grid1D:
    """
    Uniform 1D grid on an open interval (0, L), excluding boundaries.
    Dirichlet boundaries (e.g. psi=0 at 0 and L) are handled by the physics,
    not stored as degrees of freedom.
    """

    x: np.ndarray  # shape (N,), internal points only

    @property
    def dx(self) -> float:
        if self.x.size < 2:
            raise ValueError("Need at least 2 grid points to define dx.")
        return self.x[1] - self.x[0]


@dataclass
class Wavefunction1D:
    """
    Discrete wavefunction psi(x) sampled on a Grid1D.
    """

    grid: Grid1D
    values: np.ndarray  # shape (N,)

    def copy(self) -> Wavefunction1D:
        return Wavefunction1D(self.grid, self.values.copy())

    def normalize(self) -> Wavefunction1D:
        """
        Normalize with respect to the discrete L2 norm approximating
        the continuum integral ∫ |psi(x)|^2 dx on (0, L).
        """
        psi2 = np.abs(self.values) ** 2
        dx = self.grid.dx
        norm = np.sqrt(np.sum(psi2) * dx)
        if norm == 0.0:
            raise ValueError("Cannot normalize a zero wavefunction.")
        self.values /= norm
        return self

    def inner(self, other: Wavefunction1D) -> complex:
        """
        Discrete inner product ⟨self | other⟩ ≈ ∫ psi* phi dx.
        """
        if self.grid is not other.grid:
            raise ValueError("Inner product requires the same grid.")
        dx = self.grid.dx
        return np.vdot(self.values, other.values) * dx
