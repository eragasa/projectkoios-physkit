"""Deprecated one-dimensional particle-in-a-box model interfaces.

Use :mod:`projectkoios.physkit.qm.piab1d` for maintained implementations.
"""

from __future__ import annotations

from abc import ABC
from numbers import Real
from typing import ClassVar, Literal

import numpy as np

from projectkoios.physkit.constants import (
    SI,
    Constants,
    ConstantsGaussianCGS,
    ConstantsSI,
)
from projectkoios.physkit.core.boundaries import (
    AxisBoundaryConditions,
    DirichletBoundaryCondition,
)
from projectkoios.physkit.legacy.qm.models.piab import (
    Piab,
    PiabResult,
    PiabSolver,
)


class Piab1D(Piab):
    """
    Physical specification of a one-dimensional particle in a box.

    Parameters
    ----------
    x_lower:
        Lower boundary of the box.
    x_upper:
        Upper boundary of the box.
    mass:
        Particle mass in the selected unit system.
    constants:
        Physical constants expressed in the same unit system as the model
        parameters.

    Attributes
    ----------
    dimension:
        Number of spatial dimensions.
    x_lower:
        Lower boundary of the box.
    x_upper:
        Upper boundary of the box.
    length:
        Length of the box.
    mass:
        Particle mass.
    constants:
        Constants object defining the model's unit system.
    boundary_conditions:
        Homogeneous Dirichlet conditions at both boundaries.

    Notes
    -----
    The physical domain is

    .. math::

        \\Omega
        =
        [x_{\\mathrm{lower}},x_{\\mathrm{upper}}].

    The boundary conditions are

    .. math::

        \\psi(x_{\\mathrm{lower}})
        =
        \\psi(x_{\\mathrm{upper}})
        =
        0.

    The model stores only the physical specification. Numerical grids,
    operators, Hamiltonian matrices, and eigensolvers are supplied by
    solver implementations.

    All dimensional parameters must be expressed consistently in the unit
    system represented by ``constants``.
    """

    dimension: ClassVar[Literal[1]] = 1

    def __init__(
        self,
        x_lower: float,
        x_upper: float,
        mass: float,
        *,
        constants: Constants = SI,
    ) -> None:
        if isinstance(x_lower, (bool, np.bool_)) or not isinstance(x_lower, Real):
            raise TypeError("x_lower must be a real scalar.")

        if isinstance(x_upper, (bool, np.bool_)) or not isinstance(x_upper, Real):
            raise TypeError("x_upper must be a real scalar.")

        if isinstance(mass, (bool, np.bool_)) or not isinstance(mass, Real):
            raise TypeError("mass must be a real scalar.")

        if not isinstance(
            constants,
            (
                ConstantsSI,
                ConstantsGaussianCGS,
            ),
        ):
            raise TypeError(
                "constants must be an instance of ConstantsSI or ConstantsGaussianCGS."
            )

        self.x_lower: float = float(x_lower)
        self.x_upper: float = float(x_upper)
        self.mass: float = float(mass)
        self.constants: Constants = constants

        self.check_args()

        self.length: float = self.x_upper - self.x_lower

        self.boundary_conditions: AxisBoundaryConditions = AxisBoundaryConditions(
            lower=DirichletBoundaryCondition(
                value=0.0,
            ),
            upper=DirichletBoundaryCondition(
                value=0.0,
            ),
        )

    def check_args(self) -> None:
        """
        Validate the physical-model parameters.

        Raises
        ------
        ValueError
            If a coordinate, mass, or required physical constant is
            nonfinite or physically invalid.
        """
        if not np.isfinite(self.x_lower):
            raise ValueError("x_lower must be finite.")

        if not np.isfinite(self.x_upper):
            raise ValueError("x_upper must be finite.")

        if self.x_upper <= self.x_lower:
            raise ValueError("x_upper must be greater than x_lower.")

        if not np.isfinite(self.mass):
            raise ValueError("mass must be finite.")

        if self.mass <= 0.0:
            raise ValueError("mass must be greater than zero.")

        if not np.isfinite(self.constants.hbar):
            raise ValueError("constants.hbar must be finite.")

        if self.constants.hbar <= 0.0:
            raise ValueError("constants.hbar must be greater than zero.")

    def __repr__(self) -> str:
        """
        Return an unambiguous representation of the model.
        """
        return (
            f"{type(self).__name__}("
            f"x_lower={self.x_lower!r}, "
            f"x_upper={self.x_upper!r}, "
            f"mass={self.mass!r}, "
            f"constants="
            f"{type(self.constants).__name__})"
        )


class Piab1DResult(
    PiabResult,
    ABC,
):
    """
    Base result produced by a Piab1D solver.

    TISE and TDSE result classes inherit from this class.
    """

    pass


class Piab1DSolver[Piab1DResultType: Piab1DResult](
    PiabSolver[Piab1DResultType],
    ABC,
):
    """
    Base interface for a solver operating on a Piab1D model.

    Analytical TISE, finite-difference TISE, and TDSE solvers inherit from
    this class.
    """

    pass
