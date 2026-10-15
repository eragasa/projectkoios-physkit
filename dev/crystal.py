from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _as_vector(
    value: ArrayLike,
    *,
    dimension: int,
    name: str,
) -> FloatArray:
    array = np.asarray(value, dtype=np.float64)

    expected_shape = (dimension,)

    if array.shape != expected_shape:
        raise ValueError(f"{name} must have shape {expected_shape}, got {array.shape}.")

    return array


# =============================================================================
# Direct lattice
# =============================================================================


@dataclass(frozen=True)
class DirectLattice:
    """
    Base data representation of a direct lattice.

    Attributes
    ----------
    primitive_vectors:
        Array of shape (dimension, dimension).

        The rows are the primitive lattice vectors.
    """

    primitive_vectors: FloatArray

    dimension: ClassVar[int]

    def __post_init__(self) -> None:
        primitive_vectors = np.asarray(
            self.primitive_vectors,
            dtype=np.float64,
        )

        expected_shape = (self.dimension, self.dimension)

        if primitive_vectors.shape != expected_shape:
            raise ValueError(
                f"{self.__class__.__name__}.primitive_vectors must have "
                f"shape {expected_shape}, got {primitive_vectors.shape}."
            )

        object.__setattr__(self, "primitive_vectors", primitive_vectors)


@dataclass(frozen=True)
class DirectLattice1D(DirectLattice):
    """
    One-dimensional direct lattice.
    """

    dimension: ClassVar[int] = 1

    @property
    def a1(self) -> FloatArray:
        return self.primitive_vectors[0]

    @classmethod
    def from_lattice_vectors(
        cls,
        a1: ArrayLike,
    ) -> DirectLattice1D:
        a1 = _as_vector(a1, dimension=1, name="a1")
        primitive_vectors = np.vstack([a1])
        return cls(primitive_vectors=primitive_vectors)


@dataclass(frozen=True)
class DirectLattice2D(DirectLattice):
    """
    Two-dimensional direct lattice.
    """

    dimension: ClassVar[int] = 2

    @property
    def a1(self) -> FloatArray:
        return self.primitive_vectors[0]

    @property
    def a2(self) -> FloatArray:
        return self.primitive_vectors[1]

    @classmethod
    def from_lattice_vectors(
        cls,
        a1: ArrayLike,
        a2: ArrayLike,
    ) -> DirectLattice2D:
        a1 = _as_vector(a1, dimension=2, name="a1")
        a2 = _as_vector(a2, dimension=2, name="a2")

        primitive_vectors = np.vstack([a1, a2])

        return cls(primitive_vectors=primitive_vectors)


@dataclass(frozen=True)
class DirectLattice3D(DirectLattice):
    """
    Three-dimensional direct lattice.
    """

    dimension: ClassVar[int] = 3

    @property
    def a1(self) -> FloatArray:
        return self.primitive_vectors[0]

    @property
    def a2(self) -> FloatArray:
        return self.primitive_vectors[1]

    @property
    def a3(self) -> FloatArray:
        return self.primitive_vectors[2]

    @classmethod
    def from_lattice_vectors(
        cls,
        a1: ArrayLike,
        a2: ArrayLike,
        a3: ArrayLike,
    ) -> DirectLattice3D:
        a1 = _as_vector(a1, dimension=3, name="a1")
        a2 = _as_vector(a2, dimension=3, name="a2")
        a3 = _as_vector(a3, dimension=3, name="a3")

        primitive_vectors = np.vstack([a1, a2, a3])

        return cls(primitive_vectors=primitive_vectors)


# =============================================================================
# Atomic basis
# =============================================================================


@dataclass(frozen=True)
class AtomicBasis:
    """
    Base data representation of an atomic basis.

    Attributes
    ----------
    species:
        Tuple of chemical species labels.

    fractional_positions:
        Fractional coordinates of the basis atoms.

        In dimension d, this array has shape

            (n_sites, d)

        Each row stores the fractional coordinate of one basis atom.
    """

    species: tuple[str, ...]
    fractional_positions: FloatArray

    dimension: ClassVar[int]

    def __post_init__(self) -> None:
        species = tuple(self.species)

        fractional_positions = np.asarray(
            self.fractional_positions,
            dtype=np.float64,
        )

        if fractional_positions.ndim != 2:
            raise ValueError(
                f"{self.__class__.__name__}.fractional_positions must be "
                f"a 2D array, got ndim={fractional_positions.ndim}."
            )

        expected_shape_suffix = self.dimension

        if fractional_positions.shape[1] != expected_shape_suffix:
            raise ValueError(
                f"{self.__class__.__name__}.fractional_positions must have "
                f"shape (n_sites, {expected_shape_suffix}), "
                f"got {fractional_positions.shape}."
            )

        n_sites = fractional_positions.shape[0]

        if len(species) != n_sites:
            raise ValueError(
                f"{self.__class__.__name__}.species length must match "
                f"the number of basis sites. Got len(species)={len(species)} "
                f"and n_sites={n_sites}."
            )

        object.__setattr__(self, "species", species)
        object.__setattr__(self, "fractional_positions", fractional_positions)


@dataclass(frozen=True)
class AtomicBasis1D(AtomicBasis):
    """
    One-dimensional atomic basis.

    fractional_positions has shape (n_sites, 1).
    """

    dimension: ClassVar[int] = 1


@dataclass(frozen=True)
class AtomicBasis2D(AtomicBasis):
    """
    Two-dimensional atomic basis.

    fractional_positions has shape (n_sites, 2).
    """

    dimension: ClassVar[int] = 2


@dataclass(frozen=True)
class AtomicBasis3D(AtomicBasis):
    """
    Three-dimensional atomic basis.

    fractional_positions has shape (n_sites, 3).
    """

    dimension: ClassVar[int] = 3


# =============================================================================
# Crystal structure
# =============================================================================


@dataclass(frozen=True)
class CrystalStructure:
    """
    Base data representation of a crystal structure.

    A crystal structure is a direct lattice together with an atomic basis.

        CrystalStructure = DirectLattice + AtomicBasis

    This class stores only structure data. It does not generate supercells,
    reciprocal lattices, Cartesian positions, neighbor lists, or Hamiltonians.
    """

    direct_lattice: DirectLattice
    atomic_basis: AtomicBasis

    dimension: ClassVar[int]

    def __post_init__(self) -> None:
        if self.direct_lattice.dimension != self.dimension:
            raise ValueError(
                f"{self.__class__.__name__} requires a "
                f"{self.dimension}D direct lattice, got "
                f"{self.direct_lattice.dimension}D."
            )

        if self.atomic_basis.dimension != self.dimension:
            raise ValueError(
                f"{self.__class__.__name__} requires a "
                f"{self.dimension}D atomic basis, got "
                f"{self.atomic_basis.dimension}D."
            )


@dataclass(frozen=True)
class CrystalStructure1D(CrystalStructure):
    """
    One-dimensional crystal structure.
    """

    direct_lattice: DirectLattice1D
    atomic_basis: AtomicBasis1D

    dimension: ClassVar[int] = 1


@dataclass(frozen=True)
class CrystalStructure2D(CrystalStructure):
    """
    Two-dimensional crystal structure.
    """

    direct_lattice: DirectLattice2D
    atomic_basis: AtomicBasis2D

    dimension: ClassVar[int] = 2


@dataclass(frozen=True)
class CrystalStructure3D(CrystalStructure):
    """
    Three-dimensional crystal structure.
    """

    direct_lattice: DirectLattice3D
    atomic_basis: AtomicBasis3D

    dimension: ClassVar[int] = 3
