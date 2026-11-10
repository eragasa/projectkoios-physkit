r"""Define directed bonds and explicit physical geometry for ``sp3s*`` models.

A directed bond retains three conceptually separate values:

.. math::

   b=(R,s,t),\qquad R\in\mathbb Z^3,

and a composing physical geometry adds

.. math::

   \Delta r\in\mathbb R^3,\qquad
   d=\frac{\Delta r}{\lVert\Delta r\rVert_2}.

The ordinary-cell displacement :math:`R` controls coefficient placement and later
Bloch phase. The physical Cartesian displacement :math:`\Delta r` alone controls the
Slater--Koster direction. This module does not infer one from the other, from center
labels, from shell names, or from lattice metadata.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.units import (
    PhysicalUnit,
    PintUnitConverter,
    UnitSystem,
    VectorQuantity,
)


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterBondDirection:
    r"""Own one normalized three-dimensional Cartesian bond direction.

    Parameters
    ----------
    components
        Immutable nondimensional vector of shape ``(3,)`` containing
        :math:`d=(l,m,n)` from source center to target center.

    Raises
    ------
    TypeError
        If ``components`` is not exactly `VectorQuantity`.
    ValueError
        If the unit is not exactly nondimensional, the shape is not ``(3,)``, or the
        Euclidean norm differs from one by more than ``128 * eps`` in absolute value.

    Notes
    -----
    The three Cartesian components enter the represented formulas as direction
    cosines. This value owns no bond length, lattice displacement, center identity, or
    material interpretation.
    """

    components: VectorQuantity

    def __post_init__(self) -> None:
        """Validate exact nondimensional quantity ownership and unit norm."""
        if type(self.components) is not VectorQuantity:
            raise TypeError("components must be VectorQuantity")
        if self.components.unit != UnitSystem.NONDIMENSIONAL.length_unit:
            raise ValueError("components must be nondimensional")
        if self.components.magnitude.shape != (3,):
            raise ValueError("components must have shape (3,)")

        # The tolerance recognizes final binary64 normalization roundoff while still
        # rejecting approximate direction guesses or hidden magnitude information.
        tolerance = 128.0 * np.finfo(np.float64).eps
        if not math.isclose(
            float(np.linalg.norm(self.components.magnitude)),
            1.0,
            rel_tol=0.0,
            abs_tol=tolerance,
        ):
            raise ValueError("components must have unit Euclidean norm")


@dataclass(frozen=True, slots=True, kw_only=True)
class SlaterKosterDirectedBond:
    r"""Identify one directed ordinary-cell bond without Cartesian geometry.

    Parameters
    ----------
    ordinary_cell_displacement
        Exact three-dimensional integer lattice displacement :math:`R`.
    source_center_identifier
        Nonempty source-center identifier.
    target_center_identifier
        Nonempty target-center identifier.

    Raises
    ------
    TypeError
        If a field lacks its exact nominal owner.
    ValueError
        If the lattice is not three-dimensional, an identifier is empty, or the value
        would represent a zero-cell self block reserved for onsite ownership.

    Notes
    -----
    The reverse is

    .. math::

       (R,s,t)^{\mathrm{rev}}=(-R,t,s).

    The record controls ordinary-cell support only. It does not infer the physical
    vector :math:`\Delta r` or a neighbor-shell meaning.
    """

    ordinary_cell_displacement: LatticeDisplacement
    source_center_identifier: str
    target_center_identifier: str

    def __post_init__(self) -> None:
        """Validate one directed three-dimensional ordinary-cell bond."""
        if type(self.ordinary_cell_displacement) is not LatticeDisplacement:
            raise TypeError("ordinary_cell_displacement must be LatticeDisplacement")
        if self.ordinary_cell_displacement.dimension is not LatticeDimension.THREE:
            raise ValueError("sp3s* bonds require a three-dimensional lattice")

        # Endpoint identifiers remain descriptive client-owned labels. Assembly later
        # correlates them with one exact structured basis.
        for value, name in (
            (self.source_center_identifier, "source_center_identifier"),
            (self.target_center_identifier, "target_center_identifier"),
        ):
            if type(value) is not str:
                raise TypeError(f"{name} must be a built-in string")
            if not value:
                raise ValueError(f"{name} must be nonempty")

        # A zero-cell block between distinct centers remains a valid intercenter bond.
        # Only a zero-cell self block is excluded because onsite terms are independent.
        if (
            self.ordinary_cell_displacement.is_zero
            and self.source_center_identifier == self.target_center_identifier
        ):
            raise ValueError(
                "zero-cell self blocks belong to explicit onsite ownership"
            )

    def reversed(self) -> SlaterKosterDirectedBond:
        r"""Return the explicit Hermitian-reverse directed bond.

        Returns
        -------
        SlaterKosterDirectedBond
            New immutable value :math:`(-R,t,s)`.

        Raises
        ------
        RuntimeError
            If the validated displacement dimensionality is observed to have changed.
        """
        components = self.ordinary_cell_displacement.components
        if len(components) != 3:
            raise RuntimeError("validated three-dimensional displacement changed")

        # Reversal negates only the ordinary-cell displacement and swaps endpoints.
        # Cartesian geometry is deliberately owned by a separate value.
        return SlaterKosterDirectedBond(
            ordinary_cell_displacement=LatticeDisplacement(
                dimension=LatticeDimension.THREE,
                components=(-components[0], -components[1], -components[2]),
            ),
            source_center_identifier=self.target_center_identifier,
            target_center_identifier=self.source_center_identifier,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarBondGeometry:
    r"""Bind one directed bond to its physical Cartesian displacement.

    Parameters
    ----------
    directed_bond
        Exact ordinary-cell source-to-target bond identity.
    cartesian_displacement
        Immutable nonzero physical source-to-target vector :math:`\Delta r` of shape
        ``(3,)``.
    unit_system
        Physical unit-system view whose exact length unit is used by
        ``cartesian_displacement``. The default is `UnitSystem.ATOMIC`.

    Raises
    ------
    TypeError
        If values lack exact nominal owners.
    ValueError
        If the selected system is nondimensional, the vector is not a physical length,
        its unit does not exactly match the selected length unit, its shape is not
        ``(3,)``, or its norm is zero.

    Notes
    -----
    The normalized direction is derived only from

    .. math::

       d=\Delta r/\lVert\Delta r\rVert_2.

    The geometry owner remains responsible for any external crystallographic relation
    between :math:`R`, center offsets, lattice vectors, and :math:`\Delta r`.
    """

    directed_bond: SlaterKosterDirectedBond
    cartesian_displacement: VectorQuantity
    unit_system: UnitSystem = UnitSystem.ATOMIC

    def __post_init__(self) -> None:
        """Validate exact directed-bond and physical-length geometry."""
        if type(self.directed_bond) is not SlaterKosterDirectedBond:
            raise TypeError("directed_bond must be SlaterKosterDirectedBond")
        if type(self.unit_system) is not UnitSystem:
            raise TypeError("unit_system must be UnitSystem")
        if self.unit_system is UnitSystem.NONDIMENSIONAL:
            raise ValueError("unit_system must provide physical length units")
        if type(self.cartesian_displacement) is not VectorQuantity:
            raise TypeError("cartesian_displacement must be VectorQuantity")

        # Dimensional compatibility prevents an energy or reciprocal vector from being
        # normalized into a superficially valid direction.
        if type(self.cartesian_displacement.unit) is not PhysicalUnit:
            raise ValueError("cartesian_displacement must carry a physical length unit")
        if not PintUnitConverter().compatible(
            left=self.cartesian_displacement.unit,
            right=PhysicalUnit(expression="meter"),
        ):
            raise ValueError("cartesian_displacement unit must represent length")

        # Conversion is explicit and upstream. This record stores one coherent selected
        # numerical scale rather than silently converting compatible units.
        if self.cartesian_displacement.unit != self.unit_system.length_unit:
            raise ValueError(
                "cartesian_displacement must use the selected unit-system length unit"
            )
        if self.cartesian_displacement.magnitude.shape != (3,):
            raise ValueError("cartesian_displacement must have shape (3,)")
        if float(np.linalg.norm(self.cartesian_displacement.magnitude)) == 0.0:
            raise ValueError("cartesian_displacement must be nonzero")

    @property
    def bond_direction(self) -> SlaterKosterBondDirection:
        r"""Return the normalized direction derived from the physical vector.

        Returns
        -------
        SlaterKosterBondDirection
            Immutable nondimensional source-to-target direction
            :math:`\Delta r/\lVert\Delta r\rVert_2`.

        Notes
        -----
        The property performs no inference from ordinary-cell displacement or center
        embedding.
        """
        displacement = self.cartesian_displacement.magnitude

        # Normalization removes only physical length scale. The direction retains the
        # exact Cartesian orientation of the supplied vector.
        direction = displacement / np.linalg.norm(displacement)
        return SlaterKosterBondDirection(
            components=VectorQuantity(
                magnitude=np.asarray(direction, dtype=np.float64),
                unit=UnitSystem.NONDIMENSIONAL.length_unit,
            )
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarParameterizedBond:
    r"""Bind one bond geometry to an explicit parameter-sharing group.

    Parameters
    ----------
    geometry
        Validated homopolar ``sp3s*`` directed physical geometry.
    parameter_group_identifier
        Nonempty built-in string resolved by an owning model definition.

    Raises
    ------
    TypeError
        If either field lacks its exact nominal owner.
    ValueError
        If the group identifier is empty.

    Notes
    -----
    The mapping :math:`b\mapsto g(b)` expresses coefficient sharing only. It does not
    intrinsically mean a neighbor shell, material, symmetry orbit, fitted parameter
    family, or accepted scientific result.
    """

    geometry: HomopolarSp3sStarBondGeometry
    parameter_group_identifier: str

    def __post_init__(self) -> None:
        """Validate exact geometry ownership and explicit group identifier."""
        if type(self.geometry) is not HomopolarSp3sStarBondGeometry:
            raise TypeError("geometry must be HomopolarSp3sStarBondGeometry")
        if type(self.parameter_group_identifier) is not str:
            raise TypeError("parameter_group_identifier must be a built-in string")
        if not self.parameter_group_identifier:
            raise ValueError("parameter_group_identifier must be nonempty")
