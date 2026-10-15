# physkit/core/operator.py

from abc import ABC, abstractmethod

import numpy as np

from .state import Grid1D, Wavefunction1D


class LinearOperator1D(ABC):
    """
    Abstract base for linear operators acting on Wavefunction1D
    defined on a given Grid1D.
    """

    def __init__(self, grid: Grid1D):
        self.grid = grid
        self._matrix = None  # lazy matrix construction

    @abstractmethod
    def _build_matrix(self) -> np.ndarray:
        """
        Construct the matrix representation in the chosen basis.
        Must return an (N, N) ndarray.
        """
        ...

    @property
    def matrix(self) -> np.ndarray:
        if self._matrix is None:
            self._matrix = self._build_matrix()
        return self._matrix

    def apply(self, psi: Wavefunction1D) -> Wavefunction1D:
        """
        Apply the operator to a Wavefunction1D: |phi⟩ = A |psi⟩.
        """
        if psi.grid is not self.grid:
            raise ValueError("Operator and state must share the same grid.")
        values_out = self.matrix @ psi.values
        return Wavefunction1D(self.grid, values_out)
