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
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Literal, overload

from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    ModelSystemQuantity,
    ModelSystemUnit,
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


class BaseBlockHoppingTerm(ABC):
    r"""Abstract matrix-valued hopping coefficient on an ordinary lattice.

    Concrete subclasses own the coefficient representation while this class owns
    the representation-independent displacement, square-block, basis-size, unit,
    and storage-kind contract.

    Notes
    -----
    The abstraction represents one :math:`T(\mathbf R)` coefficient. It does not
    infer a conjugate partner, impose Hermiticity, or convert between dense and
    sparse storage.
    """

    __slots__ = ()

    @property
    @abstractmethod
    def displacement(self) -> LatticeDisplacement:
        """Return the exact ordinary integer-lattice displacement."""

    @property
    @abstractmethod
    def coefficient_shape(self) -> tuple[int, int]:
        """Return the coefficient's exact row and column counts."""

    @property
    def basis_size(self) -> int:
        """Return the retained primitive-cell basis dimension."""
        return self.coefficient_shape[0]

    @property
    @abstractmethod
    def energy_unit(self) -> ModelSystemUnit:
        """Return the explicit coefficient energy unit."""

    @property
    @abstractmethod
    def storage_kind(self) -> Literal["dense", "sparse"]:
        """Return the exact internal coefficient-storage representation."""

    def _validate_block_term(self) -> None:
        """Validate representation-independent block-term invariants.

        Raises
        ------
        TypeError
            If the displacement is not the exact PhysKit displacement record.
        ValueError
            If the coefficient shape is empty or nonsquare.
        """
        # An arbitrary integer tuple lacks the dimension identity required by
        # downstream model correlation, so nominal displacement ownership is exact.
        if type(self.displacement) is not LatticeDisplacement:
            raise TypeError("displacement must be LatticeDisplacement")

        # A hopping block maps one retained primitive-cell basis to that same basis.
        row_count, column_count = self.coefficient_shape
        if row_count == 0 or row_count != column_count:
            raise ValueError("coefficient must be a nonempty square matrix")


@dataclass(frozen=True, slots=True, eq=False)
class DenseBlockHoppingTerm(BaseBlockHoppingTerm):
    """Dense matrix-valued translation-invariant hopping coefficient.

    Parameters
    ----------
    displacement : LatticeDisplacement
        Ordinary primitive-cell displacement from row cell to column cell.
    coefficient : ComplexMatrixQuantity
        Immutable finite dense coefficient matrix with an explicit energy unit.

    Raises
    ------
    TypeError
        If either argument lacks its exact nominal PhysKit owner.
    ValueError
        If the coefficient is empty or nonsquare.

    Notes
    -----
    Rows and columns refer to the same ordered primitive-cell basis declared by the
    owning dense model.
    """

    displacement: LatticeDisplacement = field()
    coefficient: ComplexMatrixQuantity = field()

    def __post_init__(self) -> None:
        """Validate exact dense coefficient ownership and shared invariants."""
        if type(self.coefficient) is not ComplexMatrixQuantity:
            raise TypeError("coefficient must be ComplexMatrixQuantity")
        self._validate_block_term()

    @property
    def coefficient_shape(self) -> tuple[int, int]:
        """Return the dense coefficient shape."""
        row_count, column_count = self.coefficient.magnitude.shape
        return (int(row_count), int(column_count))

    @property
    def energy_unit(self) -> ModelSystemUnit:
        """Return the dense coefficient unit."""
        return self.coefficient.unit

    @property
    def storage_kind(self) -> Literal["dense"]:
        """Identify dense internal coefficient storage."""
        return "dense"


@dataclass(frozen=True, slots=True, eq=False)
class SparseBlockHoppingTerm(BaseBlockHoppingTerm):
    """Sparse matrix-valued translation-invariant hopping coefficient.

    Parameters
    ----------
    displacement : LatticeDisplacement
        Ordinary primitive-cell displacement from row cell to column cell.
    coefficient : ComplexSparseMatrixQuantity
        Immutable canonical-CSR coefficient matrix with an explicit energy unit.

    Raises
    ------
    TypeError
        If either argument lacks its exact nominal PhysKit owner.
    ValueError
        If the coefficient is empty or nonsquare.

    Notes
    -----
    Sparse storage describes entries inside the retained-basis block. Spatial
    support remains represented independently by the displacement inventory.
    """

    displacement: LatticeDisplacement = field()
    coefficient: ComplexSparseMatrixQuantity = field()

    def __post_init__(self) -> None:
        """Validate exact sparse coefficient ownership and shared invariants."""
        if type(self.coefficient) is not ComplexSparseMatrixQuantity:
            raise TypeError("coefficient must be ComplexSparseMatrixQuantity")
        self._validate_block_term()

    @property
    def coefficient_shape(self) -> tuple[int, int]:
        """Return the sparse coefficient shape."""
        return self.coefficient.shape

    @property
    def energy_unit(self) -> ModelSystemUnit:
        """Return the sparse coefficient unit."""
        return self.coefficient.unit

    @property
    def storage_kind(self) -> Literal["sparse"]:
        """Identify sparse internal coefficient storage."""
        return "sparse"


class BlockHoppingTerm:
    """Factory for exact dense or sparse block-hopping term records.

    Notes
    -----
    Selection uses exact nominal coefficient ownership. The factory never estimates
    numerical density or converts one representation into another.
    """

    __slots__ = ()

    @overload
    def create(
        self,
        displacement: LatticeDisplacement,
        coefficient: ComplexMatrixQuantity,
    ) -> DenseBlockHoppingTerm: ...

    @overload
    def create(
        self,
        displacement: LatticeDisplacement,
        coefficient: ComplexSparseMatrixQuantity,
    ) -> SparseBlockHoppingTerm: ...

    def create(
        self,
        displacement: LatticeDisplacement,
        coefficient: ModelSystemQuantity,
    ) -> DenseBlockHoppingTerm | SparseBlockHoppingTerm:
        """Construct the term selected by exact coefficient ownership.

        Parameters
        ----------
        displacement : LatticeDisplacement
            Ordinary integer-lattice displacement retained by the product.
        coefficient : ModelSystemQuantity
            Exact dense or sparse complex matrix quantity.

        Returns
        -------
        DenseBlockHoppingTerm or SparseBlockHoppingTerm
            Immutable term matching the supplied coefficient representation.

        Raises
        ------
        TypeError
            If ``coefficient`` is not a supported complex matrix quantity.
        ValueError
            If concrete-term validation rejects an empty or nonsquare block.
        """
        # Representation choice is nominal and never depends on numerical density.
        if type(coefficient) is ComplexMatrixQuantity:
            return DenseBlockHoppingTerm(displacement, coefficient)
        if type(coefficient) is ComplexSparseMatrixQuantity:
            return SparseBlockHoppingTerm(displacement, coefficient)
        raise TypeError(
            "coefficient must be ComplexMatrixQuantity or "
            "ComplexSparseMatrixQuantity"
        )


class BaseBlockHoppingModel(ABC):
    """Abstract owner of a homogeneous block-hopping model inventory.

    Concrete subclasses provide dense or sparse term ownership. This class owns the
    shared model identity, dimension, basis, unit, shape, displacement-order, and
    uniqueness invariants.

    Notes
    -----
    The model is an ordinary-lattice representation. It does not identify cyclic
    finite-DFT representatives, infer Hermiticity, or perform Bloch evaluation.
    """

    __slots__ = ()

    @property
    @abstractmethod
    def identifier(self) -> str:
        """Return the nonempty model identity."""

    @property
    @abstractmethod
    def dimension(self) -> LatticeDimension:
        """Return the ordinary-lattice spatial dimension."""

    @property
    @abstractmethod
    def terms(self) -> tuple[BaseBlockHoppingTerm, ...]:
        """Return the ordered homogeneous term inventory."""

    @property
    @abstractmethod
    def energy_unit(self) -> SolidStateUnit:
        """Return the exact unit shared by all coefficients."""

    @property
    @abstractmethod
    def energy_reference(self) -> str:
        """Return the common energy-zero identity."""

    @property
    @abstractmethod
    def basis_identifier(self) -> str:
        """Return the ordered primitive-cell basis identity."""

    @property
    def basis_size(self) -> int:
        """Return the internal retained-basis dimension."""
        return self.terms[0].basis_size

    @property
    @abstractmethod
    def storage_kind(self) -> Literal["dense", "sparse"]:
        """Return the exact internal coefficient-storage representation."""

    def _validate_block_model(
        self,
        owned_term_type: type[BaseBlockHoppingTerm],
    ) -> None:
        """Validate invariants shared by all concrete block models.

        Parameters
        ----------
        owned_term_type : type[BaseBlockHoppingTerm]
            Exact concrete term class admitted by the calling model.

        Raises
        ------
        TypeError
            If identities, dimension, terms, or units lack their required nominal
            owners.
        ValueError
            If identities are empty or term dimensions, shapes, units, or ordering
            are inconsistent.
        """
        # Human-readable identities are part of the scientific interpretation and
        # remain explicit rather than being inferred from coefficient contents.
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

        # Tuple order is observable and later defines deterministic Fourier and
        # serialization order, so arbitrary iterables are not accepted.
        if type(self.terms) is not tuple:
            raise TypeError("terms must be a tuple")
        if not self.terms:
            raise ValueError("terms must be nonempty")

        # A model cannot mix dense and sparse representations implicitly.
        if any(type(term) is not owned_term_type for term in self.terms):
            raise TypeError(
                f"terms must contain {owned_term_type.__name__} values"
            )
        # Every term acts on one lattice dimension and one retained state space.
        if any(
            term.displacement.dimension is not self.dimension for term in self.terms
        ):
            raise ValueError("term dimensions must match model dimension")
        coefficient_shape = self.terms[0].coefficient_shape
        if any(term.coefficient_shape != coefficient_shape for term in self.terms):
            raise ValueError("term coefficients must have one common square shape")

        # Conversion is an explicit quantity-layer action performed before model
        # construction; model units therefore match exactly.
        if not isinstance(self.energy_unit, PhysicalUnit | Unitless):
            raise TypeError("energy_unit must be PhysicalUnit or Unitless")
        if any(term.energy_unit != self.energy_unit for term in self.terms):
            raise ValueError("term coefficient units must match model energy_unit")

        # Canonical order gives deterministic identity and prevents duplicate
        # coefficients from claiming one lattice displacement.
        displacements = tuple(term.displacement.components for term in self.terms)
        if displacements != tuple(sorted(set(displacements))):
            raise ValueError("term displacements must be sorted and unique")


@dataclass(frozen=True, slots=True, eq=False)
class DenseBlockHoppingModel(BaseBlockHoppingModel):
    """Ordered ordinary-lattice Hamiltonian with dense coefficient blocks.

    Parameters
    ----------
    identifier : str
        Nonempty model identity.
    dimension : LatticeDimension
        Spatial dimension shared by every displacement.
    terms : tuple[DenseBlockHoppingTerm, ...]
        Nonempty lexicographically ordered unique displacement inventory.
    energy_unit : PhysicalUnit or Unitless
        Exact unit shared by every coefficient block.
    energy_reference : str
        Nonempty identity for the common energy zero.
    basis_identifier : str
        Nonempty identity for the ordered primitive-cell basis.

    Raises
    ------
    TypeError
        If a field lacks its required nominal owner or a term is not dense.
    ValueError
        If identities are empty or terms disagree in dimension, shape, unit, or
        canonical displacement order.

    Notes
    -----
    The model does not infer Hermiticity, lattice symmetry, truncation policy,
    cyclic finite-DFT equivalence, or scientific acceptance.
    """

    identifier: str = field()
    dimension: LatticeDimension = field()
    terms: tuple[DenseBlockHoppingTerm, ...] = field()
    energy_unit: SolidStateUnit = field()
    energy_reference: str = field()
    basis_identifier: str = field()

    def __post_init__(self) -> None:
        """Validate dense term ownership and all shared model invariants."""
        self._validate_block_model(DenseBlockHoppingTerm)

    @property
    def storage_kind(self) -> Literal["dense"]:
        """Identify dense internal coefficient storage."""
        return "dense"


@dataclass(frozen=True, slots=True, eq=False)
class SparseBlockHoppingModel(BaseBlockHoppingModel):
    """Ordered ordinary-lattice Hamiltonian with sparse coefficient blocks.

    Parameters
    ----------
    identifier : str
        Nonempty model identity.
    dimension : LatticeDimension
        Spatial dimension shared by every displacement.
    terms : tuple[SparseBlockHoppingTerm, ...]
        Nonempty lexicographically ordered unique displacement inventory.
    energy_unit : PhysicalUnit or Unitless
        Exact unit shared by every coefficient block.
    energy_reference : str
        Nonempty identity for the common energy zero.
    basis_identifier : str
        Nonempty identity for the ordered primitive-cell basis.

    Raises
    ------
    TypeError
        If a field lacks its required nominal owner or a term is not sparse.
    ValueError
        If identities are empty or terms disagree in dimension, shape, unit, or
        canonical displacement order.

    Notes
    -----
    Coefficients retain canonical immutable CSR storage. The model provides no dense
    fallback and does not infer Hermiticity or truncation policy.
    """

    identifier: str = field()
    dimension: LatticeDimension = field()
    terms: tuple[SparseBlockHoppingTerm, ...] = field()
    energy_unit: SolidStateUnit = field()
    energy_reference: str = field()
    basis_identifier: str = field()

    def __post_init__(self) -> None:
        """Validate sparse term ownership and all shared model invariants."""
        self._validate_block_model(SparseBlockHoppingTerm)

    @property
    def storage_kind(self) -> Literal["sparse"]:
        """Identify sparse internal coefficient storage."""
        return "sparse"


class BlockHoppingModel:
    """Factory for homogeneous dense or sparse block-hopping models.

    Notes
    -----
    The factory preserves supplied term order and performs no representation
    conversion. Concrete models retain the common abstract model contract.
    """

    __slots__ = ()

    @overload
    def create(
        self,
        identifier: str,
        dimension: LatticeDimension,
        terms: tuple[DenseBlockHoppingTerm, ...],
        energy_unit: SolidStateUnit,
        energy_reference: str,
        basis_identifier: str,
    ) -> DenseBlockHoppingModel: ...

    @overload
    def create(
        self,
        identifier: str,
        dimension: LatticeDimension,
        terms: tuple[SparseBlockHoppingTerm, ...],
        energy_unit: SolidStateUnit,
        energy_reference: str,
        basis_identifier: str,
    ) -> SparseBlockHoppingModel: ...

    def create(
        self,
        identifier: str,
        dimension: LatticeDimension,
        terms: tuple[DenseBlockHoppingTerm, ...]
        | tuple[SparseBlockHoppingTerm, ...],
        energy_unit: SolidStateUnit,
        energy_reference: str,
        basis_identifier: str,
    ) -> DenseBlockHoppingModel | SparseBlockHoppingModel:
        """Construct the model selected by homogeneous term ownership.

        Parameters
        ----------
        identifier : str
            Nonempty model identity retained by the concrete product.
        dimension : LatticeDimension
            Spatial dimension shared by all term displacements.
        terms : tuple
            Homogeneous exact dense-term or sparse-term inventory.
        energy_unit : PhysicalUnit or Unitless
            Exact unit shared by the coefficient blocks.
        energy_reference : str
            Identity for the common energy zero.
        basis_identifier : str
            Identity for the ordered primitive-cell basis.

        Returns
        -------
        DenseBlockHoppingModel or SparseBlockHoppingModel
            Concrete immutable model matching the term representation.

        Raises
        ------
        TypeError
            If ``terms`` is not a tuple or mixes unsupported representations.
        ValueError
            If the inventory is empty or concrete-model invariants fail.
        """
        # Empty inventories cannot determine a concrete storage family.
        if type(terms) is not tuple:
            raise TypeError("terms must be a tuple")
        if not terms:
            raise ValueError("terms must be nonempty")
        # Selection is nominal: a numerically sparse dense matrix stays dense.
        if all(type(term) is DenseBlockHoppingTerm for term in terms):
            return DenseBlockHoppingModel(
                identifier,
                dimension,
                terms,
                energy_unit,
                energy_reference,
                basis_identifier,
            )
        if all(type(term) is SparseBlockHoppingTerm for term in terms):
            return SparseBlockHoppingModel(
                identifier,
                dimension,
                terms,
                energy_unit,
                energy_reference,
                basis_identifier,
            )
        raise TypeError(
            "terms must contain only DenseBlockHoppingTerm values or only "
            "SparseBlockHoppingTerm values"
        )


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
