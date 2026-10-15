"""Deprecated misplaced one-dimensional particle-in-a-box model."""

from __future__ import annotations

from projectkoios.physkit.core.boundaries import (
    AxisBoundaryConditions,
    DirichletBoundaryCondition,
)
from projectkoios.physkit.core.constants import (
    Constants,
    ConstantsSI,
)
from projectkoios.physkit.legacy.qm.models.piab import Piab


class Piab1D(Piab):
    r"""
    One-dimensional particle-in-a-box model.

    The model describes a particle of mass :math:`m` confined to the
    interval

    .. math::

        x_{\mathrm{lower}}
        \leq x
        \leq
        x_{\mathrm{upper}}.

    The wavefunction satisfies homogeneous Dirichlet boundary conditions:

    .. math::

        \psi(x_{\mathrm{lower}}) = 0,

    .. math::

        \psi(x_{\mathrm{upper}}) = 0.

    The box length is

    .. math::

        L
        =
        x_{\mathrm{upper}}
        -
        x_{\mathrm{lower}}.

    Parameters
    ----------
    x_lower:
        Lower boundary of the box.
    x_upper:
        Upper boundary of the box.
    mass:
        Particle mass in the selected unit system.
    constants:
        Constants class defining the unit system.

    Attributes
    ----------
    x_lower:
        Lower boundary of the box.
    x_upper:
        Upper boundary of the box.
    mass:
        Particle mass.
    constants:
        Constants class defining the unit system.
    boundary_conditions:
        Homogeneous Dirichlet boundary conditions at both ends.

    Notes
    -----
    This class represents the physical model. It does not construct a
    grid, finite-difference operator, Hamiltonian matrix, or eigensolver.

    Analytical and finite-difference solutions are provided by separate
    solver classes operating on the same model instance.
    """

    def __init__(
        self,
        x_lower: float,
        x_upper: float,
        mass: float,
        *,
        constants: type[Constants] = ConstantsSI,
    ) -> None:
        self.x_lower: float = float(x_lower)
        self.x_upper: float = float(x_upper)
        self.mass: float = float(mass)
        self.constants: type[Constants] = constants

        self.check_args()

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
        Validate the physical model parameters.

        Raises
        ------
        TypeError
            If ``constants`` is not a constants class.
        ValueError
            If the interval has nonpositive length or the particle mass
            is not positive.
        """
        if self.x_upper <= self.x_lower:
            raise ValueError("x_upper must be greater than x_lower.")

        if self.mass <= 0.0:
            raise ValueError("mass must be positive.")

        if not isinstance(self.constants, type):
            raise TypeError("constants must be a constants class.")

        if not issubclass(
            self.constants,
            Constants,
        ):
            raise TypeError("constants must inherit from Constants.")

    @property
    def length(self) -> float:
        """
        Return the length of the box.

        Returns
        -------
        float
            Box length in the selected unit system.
        """
        return self.x_upper - self.x_lower

    @property
    def dimension(self) -> int:
        """
        Return the spatial dimension of the model.

        Returns
        -------
        int
            Spatial dimension, equal to one.
        """
        return 1

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}("
            f"x_lower={self.x_lower!r}, "
            f"x_upper={self.x_upper!r}, "
            f"mass={self.mass!r}, "
            f"constants={self.constants.__name__})"
        )
