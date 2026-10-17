"""Immutable calculator-neutral periodic unit-cell models."""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from projectkoios.physkit.periodic.lattice.lattice3d import DirectLattice3D
from projectkoios.physkit.units.quantities import (
    MatrixQuantity,
    ModelSystemUnitConverter,
    PhysicalUnit,
    ScalarQuantity,
    Unitless,
    VectorQuantity,
)

_ELEMENT_SYMBOL = re.compile(r"[A-Z][a-z]?")


@dataclass(frozen=True, slots=True)
class Atom:
    """Represent one chemical symbol at fractional unit-cell coordinates."""

    symbol: str
    position_fractional: VectorQuantity

    def __post_init__(self) -> None:
        if type(self.symbol) is not str or not _ELEMENT_SYMBOL.fullmatch(self.symbol):
            raise ValueError("atom symbol must be an element symbol")
        if type(self.position_fractional) is not VectorQuantity:
            raise TypeError("atom position_fractional must be a VectorQuantity")
        if type(self.position_fractional.unit) is not Unitless:
            raise ValueError("atom position_fractional must be explicitly unitless")
        if self.position_fractional.magnitude.shape != (3,):
            raise ValueError("atom position_fractional must contain three coordinates")


@dataclass(frozen=True, slots=True)
class AtomicBasis:
    """Represent an ordered, nonempty tuple of atoms in one unit cell."""

    atoms: tuple[Atom, ...]

    def __post_init__(self) -> None:
        if type(self.atoms) is not tuple:
            raise TypeError("atomic basis atoms must be a tuple")
        if not self.atoms:
            raise ValueError("atomic basis must contain at least one atom")
        if any(type(atom) is not Atom for atom in self.atoms):
            raise TypeError("atomic basis must contain Atom values")


@dataclass(frozen=True, slots=True)
class UnitCell:
    r"""Compose dimensionless and physical column-basis representations.

    The direct-lattice matrix is

    .. math::

        A = [\mathbf a_1\ \mathbf a_2\ \mathbf a_3].

    For positive lattice parameter ``a``, the physical cell matrix is

    .. math::

        H = aA = [\mathbf h_1\ \mathbf h_2\ \mathbf h_3].

    The supplied direct lattice is copied so that the unit cell owns immutable
    geometry independent of later caller mutation.
    """

    direct_lattice: DirectLattice3D
    lattice_parameter: ScalarQuantity
    atomic_basis: AtomicBasis
    _A: MatrixQuantity = field(init=False, repr=False)
    _H: MatrixQuantity = field(init=False, repr=False)
    _a_vectors: tuple[VectorQuantity, VectorQuantity, VectorQuantity] = field(
        init=False,
        repr=False,
    )
    _h_vectors: tuple[VectorQuantity, VectorQuantity, VectorQuantity] = field(
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.direct_lattice, DirectLattice3D):
            raise TypeError("direct_lattice must be a DirectLattice3D")
        if type(self.lattice_parameter) is not ScalarQuantity:
            raise TypeError("lattice_parameter must be a ScalarQuantity")
        if type(self.lattice_parameter.unit) is not PhysicalUnit:
            raise ValueError("lattice_parameter must have an explicit physical unit")
        if not ModelSystemUnitConverter.compatible(
            self.lattice_parameter.unit,
            PhysicalUnit("meter"),
        ):
            raise ValueError(
                "lattice_parameter must have physical length dimensionality"
            )
        if self.lattice_parameter.magnitude <= 0.0:
            raise ValueError("lattice_parameter must be positive")
        if type(self.atomic_basis) is not AtomicBasis:
            raise TypeError("atomic_basis must be an AtomicBasis")

        source_a = np.asarray(self.direct_lattice.A, dtype=np.float64)
        canonical_lattice = DirectLattice3D(
            a1=source_a[:, 0],
            a2=source_a[:, 1],
            a3=source_a[:, 2],
        )
        a_matrix = MatrixQuantity(magnitude=source_a, unit=Unitless())
        h_matrix = MatrixQuantity(
            magnitude=a_matrix.magnitude * self.lattice_parameter.magnitude,
            unit=self.lattice_parameter.unit,
        )
        a_vectors = tuple(
            VectorQuantity(
                magnitude=a_matrix.magnitude[:, index],
                unit=Unitless(),
            )
            for index in range(3)
        )
        h_vectors = tuple(
            VectorQuantity(
                magnitude=h_matrix.magnitude[:, index],
                unit=self.lattice_parameter.unit,
            )
            for index in range(3)
        )
        object.__setattr__(self, "direct_lattice", canonical_lattice)
        object.__setattr__(self, "_A", a_matrix)
        object.__setattr__(self, "_H", h_matrix)
        object.__setattr__(self, "_a_vectors", a_vectors)
        object.__setattr__(self, "_h_vectors", h_vectors)

    @property
    def A(self) -> MatrixQuantity:
        """Return the immutable dimensionless direct-lattice matrix."""
        return self._A

    @property
    def H(self) -> MatrixQuantity:
        """Return the immutable physical cell matrix."""
        return self._H

    @property
    def a1(self) -> VectorQuantity:
        """Return the first dimensionless direct-lattice vector."""
        return self._a_vectors[0]

    @property
    def a2(self) -> VectorQuantity:
        """Return the second dimensionless direct-lattice vector."""
        return self._a_vectors[1]

    @property
    def a3(self) -> VectorQuantity:
        """Return the third dimensionless direct-lattice vector."""
        return self._a_vectors[2]

    @property
    def h1(self) -> VectorQuantity:
        """Return the first physical cell vector."""
        return self._h_vectors[0]

    @property
    def h2(self) -> VectorQuantity:
        """Return the second physical cell vector."""
        return self._h_vectors[1]

    @property
    def h3(self) -> VectorQuantity:
        """Return the third physical cell vector."""
        return self._h_vectors[2]


@dataclass(frozen=True, slots=True)
class ConventionalUnitCell(UnitCell):
    """Identify a unit cell containing a declared conventional representation."""


@dataclass(frozen=True, slots=True)
class PrimitiveUnitCell(UnitCell):
    """Identify a unit cell containing a declared primitive representation."""
