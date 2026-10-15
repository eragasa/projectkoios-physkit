"""Boundary conditions for structured numerical grids."""

from __future__ import annotations

from enum import Enum
from numbers import Complex

import numpy as np

__all__ = [
    "BoundaryConditionType",
    "BoundaryCondition",
    "AxisBoundaryConditions",
    "DirichletBoundaryCondition",
    "NeumannBoundaryCondition",
    "RobinBoundaryCondition",
    "PeriodicBoundaryCondition",
    "BlochBoundaryCondition",
    "AsymptoticDecayBoundaryCondition",
    "SommerfeldRadiationBoundaryCondition",
]


class BoundaryConditionType(Enum):
    """Boundary-condition type.

    Attributes
    ----------
    DIRICHLET
        Prescribed field value.
    NEUMANN
        Prescribed outward-normal derivative.
    ROBIN
        Linear relation between field value and normal derivative.
    PERIODIC
        Periodic coupling between paired boundaries.
    BLOCH
        Bloch coupling between paired boundaries.
    ASYMPTOTIC_DECAY
        Decay condition at spatial infinity.
    SOMMERFELD_RADIATION
        Local outgoing-wave condition.
    """

    DIRICHLET = "dirichlet"
    NEUMANN = "neumann"
    ROBIN = "robin"
    PERIODIC = "periodic"
    BLOCH = "bloch"
    ASYMPTOTIC_DECAY = "asymptotic_decay"
    SOMMERFELD_RADIATION = "sommerfeld_radiation"


class BoundaryCondition:
    """Base class for boundary-condition specifications.

    Attributes
    ----------
    type : BoundaryConditionType
        Boundary-condition type.

    Notes
    -----
    A boundary-condition object stores a mathematical constraint. It does
    not identify a boundary region or construct a numerical operator.
    """

    type: BoundaryConditionType


class AxisBoundaryConditions:
    """Boundary conditions associated with one coordinate axis.

    Parameters
    ----------
    lower : BoundaryCondition
        Condition on the lower side of the coordinate interval.
    upper : BoundaryCondition
        Condition on the upper side of the coordinate interval.

    Attributes
    ----------
    lower : BoundaryCondition
        Lower-side boundary condition.
    upper : BoundaryCondition
        Upper-side boundary condition.

    Raises
    ------
    TypeError
        If either argument is not a ``BoundaryCondition``.
    ValueError
        If a paired condition occurs on only one side or if the lower and
        upper paired-condition types do not match.

    Notes
    -----
    This class assigns conditions to the two boundary faces normal to one
    axis of a structured tensor-product grid. A grid in :math:`d`
    dimensions uses one ``AxisBoundaryConditions`` object per axis.
    """

    def __init__(
        self,
        lower: BoundaryCondition,
        upper: BoundaryCondition,
    ) -> None:
        if not isinstance(lower, BoundaryCondition):
            raise TypeError("lower must be a BoundaryCondition.")

        if not isinstance(upper, BoundaryCondition):
            raise TypeError("upper must be a BoundaryCondition.")

        self.lower: BoundaryCondition = lower
        self.upper: BoundaryCondition = upper

        self.check_args()

    def check_args(self) -> None:
        """Validate the lower and upper condition pairing.

        Raises
        ------
        ValueError
            If a periodic or Bloch condition is not paired with the same
            condition type on the opposite side.
        """
        paired_types = {
            BoundaryConditionType.PERIODIC,
            BoundaryConditionType.BLOCH,
        }

        lower_is_paired = self.lower.type in paired_types
        upper_is_paired = self.upper.type in paired_types

        if lower_is_paired != upper_is_paired:
            raise ValueError(
                "Periodic and Bloch conditions must occur on both sides of an axis."
            )

        if lower_is_paired and self.lower.type is not self.upper.type:
            raise ValueError("Paired boundary-condition types must match.")

    def __repr__(self) -> str:
        return f"AxisBoundaryConditions(lower={self.lower!r}, upper={self.upper!r})"


class DirichletBoundaryCondition(BoundaryCondition):
    r"""Prescribed field value on a boundary.

    Parameters
    ----------
    value : numbers.Complex, optional
        Constant field value. The default is ``0.0``.

    Attributes
    ----------
    value : complex
        Prescribed boundary value.
    type : BoundaryConditionType
        Boundary-condition type.

    Raises
    ------
    TypeError
        If ``value`` is not a real or complex scalar.
    ValueError
        If ``value`` is not finite.

    Notes
    -----
    The condition is

    .. math::

        u\vert_{\Gamma}=g.
    """

    def __init__(
        self,
        value: Complex = 0.0,
    ) -> None:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Complex):
            raise TypeError("value must be a real or complex scalar.")

        self.value: complex = complex(value)
        self.type: BoundaryConditionType = BoundaryConditionType.DIRICHLET

        self.check_args()

    def check_args(self) -> None:
        """Validate the prescribed value."""
        if not np.isfinite(self.value):
            raise ValueError("value must be finite.")

    def __repr__(self) -> str:
        return f"DirichletBoundaryCondition(value={self.value!r})"


class NeumannBoundaryCondition(BoundaryCondition):
    r"""Prescribed outward-normal derivative on a boundary.

    Parameters
    ----------
    value : numbers.Complex, optional
        Constant outward-normal derivative. The default is ``0.0``.

    Attributes
    ----------
    value : complex
        Prescribed outward-normal derivative.
    type : BoundaryConditionType
        Boundary-condition type.

    Raises
    ------
    TypeError
        If ``value`` is not a real or complex scalar.
    ValueError
        If ``value`` is not finite.

    Notes
    -----
    The condition is

    .. math::

        \frac{\partial u}{\partial n}\bigg\vert_{\Gamma}=q,

    where :math:`\hat{\mathbf n}` is the outward unit normal.
    """

    def __init__(
        self,
        value: Complex = 0.0,
    ) -> None:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Complex):
            raise TypeError("value must be a real or complex scalar.")

        self.value: complex = complex(value)
        self.type: BoundaryConditionType = BoundaryConditionType.NEUMANN

        self.check_args()

    def check_args(self) -> None:
        """Validate the prescribed derivative."""
        if not np.isfinite(self.value):
            raise ValueError("value must be finite.")

    def __repr__(self) -> str:
        return f"NeumannBoundaryCondition(value={self.value!r})"


class RobinBoundaryCondition(BoundaryCondition):
    r"""Linear relation between field value and normal derivative.

    Parameters
    ----------
    alpha : numbers.Complex
        Coefficient multiplying the field value.
    beta : numbers.Complex
        Coefficient multiplying the outward-normal derivative.
    value : numbers.Complex, optional
        Right-hand-side value. The default is ``0.0``.

    Attributes
    ----------
    alpha : complex
        Field-value coefficient.
    beta : complex
        Normal-derivative coefficient.
    value : complex
        Right-hand-side value.
    type : BoundaryConditionType
        Boundary-condition type.

    Raises
    ------
    TypeError
        If a parameter is not a real or complex scalar.
    ValueError
        If a parameter is nonfinite or ``alpha`` and ``beta`` are both zero.

    Notes
    -----
    The condition is

    .. math::

        \alpha u
        +
        \beta\frac{\partial u}{\partial n}
        =g.
    """

    def __init__(
        self,
        alpha: Complex,
        beta: Complex,
        value: Complex = 0.0,
    ) -> None:
        if isinstance(alpha, (bool, np.bool_)) or not isinstance(alpha, Complex):
            raise TypeError("alpha must be a real or complex scalar.")

        if isinstance(beta, (bool, np.bool_)) or not isinstance(beta, Complex):
            raise TypeError("beta must be a real or complex scalar.")

        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Complex):
            raise TypeError("value must be a real or complex scalar.")

        self.alpha: complex = complex(alpha)
        self.beta: complex = complex(beta)
        self.value: complex = complex(value)
        self.type: BoundaryConditionType = BoundaryConditionType.ROBIN

        self.check_args()

    def check_args(self) -> None:
        """Validate the Robin coefficients and value."""
        if not np.isfinite(self.alpha):
            raise ValueError("alpha must be finite.")

        if not np.isfinite(self.beta):
            raise ValueError("beta must be finite.")

        if not np.isfinite(self.value):
            raise ValueError("value must be finite.")

        if self.alpha == 0.0 and self.beta == 0.0:
            raise ValueError("alpha and beta must not both be zero.")

    def __repr__(self) -> str:
        return (
            "RobinBoundaryCondition("
            f"alpha={self.alpha!r}, "
            f"beta={self.beta!r}, "
            f"value={self.value!r})"
        )


class PeriodicBoundaryCondition(BoundaryCondition):
    r"""Periodic coupling between paired boundaries.

    Attributes
    ----------
    type : BoundaryConditionType
        Boundary-condition type.

    Notes
    -----
    The paired-boundary relation is

    .. math::

        u(\mathbf r+\mathbf L)=u(\mathbf r).
    """

    def __init__(self) -> None:
        self.type: BoundaryConditionType = BoundaryConditionType.PERIODIC

    def __repr__(self) -> str:
        return "PeriodicBoundaryCondition()"


class BlochBoundaryCondition(BoundaryCondition):
    r"""Bloch coupling between paired boundaries.

    Attributes
    ----------
    type : BoundaryConditionType
        Boundary-condition type.

    Notes
    -----
    The paired-boundary relation is

    .. math::

        u_{\mathbf k}(\mathbf r+\mathbf L)
        =
        e^{i\mathbf k\cdot\mathbf L}
        u_{\mathbf k}(\mathbf r).

    This object declares the relation but does not store a fixed wavevector
    or phase. The Bloch solver supplies :math:`\mathbf k`, and the boundary
    pairing supplies :math:`\mathbf L`.
    """

    def __init__(self) -> None:
        self.type: BoundaryConditionType = BoundaryConditionType.BLOCH

    def __repr__(self) -> str:
        return "BlochBoundaryCondition()"


class AsymptoticDecayBoundaryCondition(BoundaryCondition):
    r"""Decay condition at spatial infinity.

    Attributes
    ----------
    type : BoundaryConditionType
        Boundary-condition type.

    Notes
    -----
    The asymptotic condition is

    .. math::

        \lim_{\lVert\mathbf r\rVert\rightarrow\infty}
        u(\mathbf r)=0.

    A finite-domain solver must select and validate a numerical truncation
    that approximates this condition.
    """

    def __init__(self) -> None:
        self.type: BoundaryConditionType = BoundaryConditionType.ASYMPTOTIC_DECAY

    def __repr__(self) -> str:
        return "AsymptoticDecayBoundaryCondition()"


class SommerfeldRadiationBoundaryCondition(BoundaryCondition):
    r"""Local outgoing-wave boundary condition.

    Parameters
    ----------
    wavenumber : numbers.Complex
        Wavenumber of the outgoing field.

    Attributes
    ----------
    wavenumber : complex
        Outgoing-field wavenumber.
    type : BoundaryConditionType
        Boundary-condition type.

    Raises
    ------
    TypeError
        If ``wavenumber`` is not a real or complex scalar.
    ValueError
        If ``wavenumber`` is not finite.

    Notes
    -----
    The first-order radiation condition is

    .. math::

        \frac{\partial u}{\partial n}=iku.
    """

    def __init__(
        self,
        wavenumber: Complex,
    ) -> None:
        if isinstance(wavenumber, (bool, np.bool_)) or not isinstance(
            wavenumber, Complex
        ):
            raise TypeError("wavenumber must be a real or complex scalar.")

        self.wavenumber: complex = complex(wavenumber)

        self.type: BoundaryConditionType = BoundaryConditionType.SOMMERFELD_RADIATION

        self.check_args()

    def check_args(self) -> None:
        """Validate the outgoing-field wavenumber."""
        if not np.isfinite(self.wavenumber):
            raise ValueError("wavenumber must be finite.")

    def __repr__(self) -> str:
        return f"SommerfeldRadiationBoundaryCondition(wavenumber={self.wavenumber!r})"
