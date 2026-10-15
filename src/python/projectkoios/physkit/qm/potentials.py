from abc import ABC, abstractmethod

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix


class Potential1D(ABC):
    """
    Base class for local scalar 1D potentials.

    A Potential1D represents a function V(x). In the discrete position basis,
    the corresponding potential-energy operator is diagonal:

        V_ij = V(x_i) delta_ij

    Subclasses only need to define eval().
    """

    @abstractmethod
    def eval(self, x: np.ndarray) -> np.ndarray:
        """
        Evaluate the potential on a 1D grid.

        Parameters
        ----------
        x:
            Real-space grid points.

        Returns
        -------
        np.ndarray
            Array of potential values V(x_i).
        """
        raise NotImplementedError

    def matrix(self, x: np.ndarray) -> np.ndarray:
        """
        Dense matrix representation of the potential operator.
        """
        return np.diag(self.eval(x))

    def sparse_matrix(self, x: np.ndarray) -> csr_matrix:
        """
        Sparse matrix representation of the potential operator.
        """
        return sp.diags(self.eval(x), offsets=0, format="csr")


class ConstantPotential1D(Potential1D):
    """
    Constant potential:

        V(x) = V0
    """

    def __init__(self, value: float = 0.0):
        self.value = float(value)

    def eval(self, x: np.ndarray) -> np.ndarray:
        return self.value * np.ones_like(x, dtype=float)
