"""Immutable scalar and block hopping models and localized perturbations.

These records describe reduced lattice Hamiltonian data independently of atomic crystal
structures, DFT provenance, model-class fitting, or scientific acceptance. Energy units
use PhysKit's Pint-backed unit records, and every model retains explicit basis and
energy-reference identities.

This implementation preserves the represented behavior of
``ksdft2effmass.solid_state.lattice_models`` at donor revision
``7bd913151f7e61ed2bdba593df920be36573b502`` while using PhysKit's native unit
records.

References
----------
.. [1] J. C. Slater and G. F. Koster, "Simplified LCAO Method for the
   Periodic Potential Problem," *Physical Review* 94, 1498 (1954).
   https://doi.org/10.1103/PhysRev.94.1498
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    PhysicalUnit,
    Unitless,
)

from .finite_domain import LatticeCoordinate, LatticeDimension, LatticeDisplacement

type SolidStateUnit = PhysicalUnit | Unitless
"""Explicit physical or dimensionless unit admitted by solid-state records."""


@dataclass(frozen=True, slots=True)
class ScalarHoppingTerm:
    """Represent one scalar translation-invariant hopping coefficient.

    Parameters
    ----------
    displacement
        Integer primitive-cell displacement from row cell to column cell.
    real
        Finite built-in-float real part of the coefficient.
    imaginary
        Finite built-in-float imaginary part of the coefficient.
    """

    displacement: LatticeDisplacement
    real: float
    imaginary: float

    def __post_init__(self) -> None:
        """Validate displacement and finite binary64 components."""
        if type(self.displacement) is not LatticeDisplacement:
            raise TypeError("displacement must be LatticeDisplacement")
        if type(self.real) is not float or type(self.imaginary) is not float:
            raise TypeError("hopping components must be built-in floats")
        if not math.isfinite(self.real) or not math.isfinite(self.imaginary):
            raise ValueError("hopping components must be finite")

    @property
    def value(self) -> complex:
        """Return the represented built-in complex coefficient."""
        return complex(self.real, self.imaginary)


@dataclass(frozen=True, slots=True)
class ScalarHoppingModel:
    """Represent one ordered scalar translation-invariant lattice Hamiltonian.

    Parameters
    ----------
    identifier
        Nonempty model identity.
    dimension
        Exact supported spatial dimension.
    terms
        Nonempty hopping terms sorted lexicographically by displacement with no
        duplicates. Hermiticity is not inferred by construction.
    energy_unit
        Explicit Pint-backed physical or dimensionless unit shared by every term.
    energy_reference
        Nonempty identity for the common energy zero.
    basis_identifier
        Nonempty identity for the scalar primitive-cell basis.
    """

    identifier: str
    dimension: LatticeDimension
    terms: tuple[ScalarHoppingTerm, ...]
    energy_unit: SolidStateUnit
    energy_reference: str
    basis_identifier: str

    def __post_init__(self) -> None:
        """Validate intrinsic identity, dimensional, ordering, and unit invariants."""
        for value, name in (
            (self.identifier, "identifier"),
            (self.energy_reference, "energy_reference"),
            (self.basis_identifier, "basis_identifier"),
        ):
            if type(value) is not str:
                raise TypeError(f"{name} must be a string")
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if type(self.dimension) is not LatticeDimension:
            raise TypeError("dimension must be LatticeDimension")
        if type(self.terms) is not tuple:
            raise TypeError("terms must be a tuple")
        if not self.terms:
            raise ValueError("terms must be nonempty")
        if any(type(term) is not ScalarHoppingTerm for term in self.terms):
            raise TypeError("terms must contain ScalarHoppingTerm values")
        if any(
            term.displacement.dimension is not self.dimension for term in self.terms
        ):
            raise ValueError("term dimensions must match model dimension")
        displacements = tuple(term.displacement.components for term in self.terms)
        if displacements != tuple(sorted(set(displacements))):
            raise ValueError("term displacements must be sorted and unique")
        if not isinstance(self.energy_unit, PhysicalUnit | Unitless):
            raise TypeError("energy_unit must be PhysicalUnit or Unitless")


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingTerm:
    """Represent one matrix-valued translation-invariant hopping coefficient.

    ``coefficient`` carries the energy unit explicitly. Its rows and columns refer to
    the same ordered primitive-cell basis declared by the owning model.
    """

    displacement: LatticeDisplacement
    coefficient: ComplexMatrixQuantity

    def __post_init__(self) -> None:
        """Validate one finite immutable nonempty square coefficient block."""
        if type(self.displacement) is not LatticeDisplacement:
            raise TypeError("displacement must be LatticeDisplacement")
        if type(self.coefficient) is not ComplexMatrixQuantity:
            raise TypeError("coefficient must be ComplexMatrixQuantity")
        row_count, column_count = self.coefficient.magnitude.shape
        if row_count == 0 or row_count != column_count:
            raise ValueError("coefficient must be a nonempty square matrix")


@dataclass(frozen=True, slots=True, eq=False)
class BlockHoppingModel:
    """Represent one ordered matrix-valued lattice Hamiltonian.

    The model owns displacement and coefficient order, lattice dimension, energy unit,
    energy-reference identity, and ordered primitive-cell basis identity. It does not
    infer Hermiticity, lattice symmetry, truncation policy, or scientific acceptance.
    Displacements are ordinary integer-lattice translations; no cyclic representative
    modulus or finite-DFT alias equivalence is implied.
    """

    identifier: str
    dimension: LatticeDimension
    terms: tuple[BlockHoppingTerm, ...]
    energy_unit: SolidStateUnit
    energy_reference: str
    basis_identifier: str

    def __post_init__(self) -> None:
        """Validate identity, dimensions, common matrix shape, units, and order."""
        for value, name in (
            (self.identifier, "identifier"),
            (self.energy_reference, "energy_reference"),
            (self.basis_identifier, "basis_identifier"),
        ):
            if type(value) is not str:
                raise TypeError(f"{name} must be a string")
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if type(self.dimension) is not LatticeDimension:
            raise TypeError("dimension must be LatticeDimension")
        if type(self.terms) is not tuple:
            raise TypeError("terms must be a tuple")
        if not self.terms:
            raise ValueError("terms must be nonempty")
        if any(type(term) is not BlockHoppingTerm for term in self.terms):
            raise TypeError("terms must contain BlockHoppingTerm values")
        if any(
            term.displacement.dimension is not self.dimension for term in self.terms
        ):
            raise ValueError("term dimensions must match model dimension")
        coefficient_shape = self.terms[0].coefficient.magnitude.shape
        if any(
            term.coefficient.magnitude.shape != coefficient_shape for term in self.terms
        ):
            raise ValueError("term coefficients must have one common square shape")
        if not isinstance(self.energy_unit, PhysicalUnit | Unitless):
            raise TypeError("energy_unit must be PhysicalUnit or Unitless")
        if any(term.coefficient.unit != self.energy_unit for term in self.terms):
            raise ValueError("term coefficient units must match model energy_unit")
        displacements = tuple(term.displacement.components for term in self.terms)
        if displacements != tuple(sorted(set(displacements))):
            raise ValueError("term displacements must be sorted and unique")

    @property
    def basis_size(self) -> int:
        """Return the row and column count of each coefficient block."""
        return int(self.terms[0].coefficient.magnitude.shape[0])


@dataclass(frozen=True, slots=True)
class LocalizedOnsiteTerm:
    """Represent one scalar localized onsite coefficient.

    Parameters
    ----------
    site
        Integer lattice site relative to an explicitly identified defect origin.
    real
        Finite built-in-float real part.
    imaginary
        Finite built-in-float imaginary part. Hermitian onsite perturbations require
        this value to be zero, but this structural record does not impose that policy.
    """

    site: LatticeCoordinate
    real: float
    imaginary: float

    def __post_init__(self) -> None:
        """Validate site and finite coefficient components."""
        if type(self.site) is not LatticeCoordinate:
            raise TypeError("site must be LatticeCoordinate")
        if type(self.real) is not float or type(self.imaginary) is not float:
            raise TypeError("onsite components must be built-in floats")
        if not math.isfinite(self.real) or not math.isfinite(self.imaginary):
            raise ValueError("onsite components must be finite")

    @property
    def value(self) -> complex:
        """Return the represented built-in complex coefficient."""
        return complex(self.real, self.imaginary)


@dataclass(frozen=True, slots=True)
class LocalizedBondTerm:
    """Represent one directed localized scalar bond coefficient.

    Parameters
    ----------
    start
        Integer start site relative to the defect origin.
    displacement
        Nonzero integer displacement from start to target.
    real
        Finite built-in-float real part.
    imaginary
        Finite built-in-float imaginary part.
    """

    start: LatticeCoordinate
    displacement: LatticeDisplacement
    real: float
    imaginary: float

    def __post_init__(self) -> None:
        """Validate dimensions, nonzero displacement, and finite coefficient."""
        if type(self.start) is not LatticeCoordinate:
            raise TypeError("start must be LatticeCoordinate")
        if type(self.displacement) is not LatticeDisplacement:
            raise TypeError("displacement must be LatticeDisplacement")
        if self.start.dimension is not self.displacement.dimension:
            raise ValueError("start and displacement dimensions must agree")
        if self.displacement.is_zero:
            raise ValueError("bond displacement must be nonzero")
        if type(self.real) is not float or type(self.imaginary) is not float:
            raise TypeError("bond components must be built-in floats")
        if not math.isfinite(self.real) or not math.isfinite(self.imaginary):
            raise ValueError("bond components must be finite")

    @property
    def value(self) -> complex:
        """Return the represented built-in complex coefficient."""
        return complex(self.real, self.imaginary)


type LocalizedTerm = LocalizedOnsiteTerm | LocalizedBondTerm
"""Closed scalar localized-term union."""


@dataclass(frozen=True, slots=True)
class LocalizedPerturbation:
    """Represent one deterministic finite localized scalar perturbation.

    Parameters
    ----------
    identifier
        Nonempty perturbation identity.
    dimension
        Exact supported spatial dimension.
    terms
        Nonempty tuple of onsite and directed bond terms. Terms must be in canonical
        order: onsite terms first by site, then bond terms by start and displacement.
        Duplicate represented terms are rejected.
    energy_unit
        Explicit Pint-backed unit shared by every coefficient.
    energy_reference
        Nonempty identity for the common energy zero.
    basis_identifier
        Nonempty scalar basis identity.
    """

    identifier: str
    dimension: LatticeDimension
    terms: tuple[LocalizedTerm, ...]
    energy_unit: SolidStateUnit
    energy_reference: str
    basis_identifier: str

    def __post_init__(self) -> None:
        """Validate identity, dimensional compatibility, canonical order, and unit."""
        for value, name in (
            (self.identifier, "identifier"),
            (self.energy_reference, "energy_reference"),
            (self.basis_identifier, "basis_identifier"),
        ):
            if type(value) is not str:
                raise TypeError(f"{name} must be a string")
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if type(self.dimension) is not LatticeDimension:
            raise TypeError("dimension must be LatticeDimension")
        if type(self.terms) is not tuple:
            raise TypeError("terms must be a tuple")
        if not self.terms:
            raise ValueError("terms must be nonempty")
        keys: list[tuple[int, ...]] = []
        for term in self.terms:
            if type(term) is LocalizedOnsiteTerm:
                if term.site.dimension is not self.dimension:
                    raise ValueError("onsite term dimension must match perturbation")
                keys.append((0, *term.site.components))
            elif type(term) is LocalizedBondTerm:
                if term.start.dimension is not self.dimension:
                    raise ValueError("bond term dimension must match perturbation")
                keys.append((1, *term.start.components, *term.displacement.components))
            else:
                raise TypeError("terms must contain supported localized term values")
        if tuple(keys) != tuple(sorted(set(keys))):
            raise ValueError("localized terms must be canonically sorted and unique")
        if not isinstance(self.energy_unit, PhysicalUnit | Unitless):
            raise TypeError("energy_unit must be PhysicalUnit or Unitless")
