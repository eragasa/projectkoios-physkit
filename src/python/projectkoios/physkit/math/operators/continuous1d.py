# continuous1d.py
# ------------------------------------------------------------
# 1D Continuous Linear Operators for PDEs
#
# This module defines abstract and concrete differential
# operators acting on scalar fields phi(x) represented as
# SymPy Functions. No discretization is performed here.
#
# Namespace recommendation:
#    math.operators.continuous.continuous1d
# ------------------------------------------------------------

from abc import ABC, abstractmethod

import sympy as sp

_DEFAULT_X = sp.Symbol("x")


# ------------------------------------------------------------
# Base continuous operator
# ------------------------------------------------------------
class ContinuousOperator1D(ABC):
    """
    Abstract base class for continuous 1D linear differential operators.

    Operators act on symbolic functions phi(x) via:
        result(x) = self.apply(phi)

    where phi is a SymPy Function, not a discrete array.
    """

    def __init__(self, x_symbol=_DEFAULT_X):
        self.x = x_symbol

    @abstractmethod
    def apply(self, phi):
        """
        Apply the operator to a SymPy function phi(x).

        Parameters
        ----------
        phi : sympy.Function
            Scalar field phi(x) on which the operator acts.

        Returns
        -------
        sympy expression
            The symbolic expression representing the operator applied to phi.
        """
        raise NotImplementedError


# ------------------------------------------------------------
# Poisson operator P[phi] = d/dx ( eps(x) * dphi/dx )
# ------------------------------------------------------------
class PoissonOperator1D(ContinuousOperator1D):
    """
    Continuous Poisson operator:

        P[phi] = d/dx ( eps(x) * dphi/dx )

    where eps(x) is a SymPy expression.

    This operator does not include boundary conditions or
    discretization; it is purely symbolic.
    """

    def __init__(self, eps_expr, x_symbol=_DEFAULT_X):
        super().__init__(x_symbol=x_symbol)

        # eps(x) must be a SymPy expression depending on x
        self.eps = eps_expr

    def apply(self, phi):
        """
        Apply the Poisson operator to phi(x):

            return d/dx [ eps(x) * dphi/dx ].
        """
        x = self.x
        eps = self.eps
        return sp.diff(eps * sp.diff(phi(x), x), x)
