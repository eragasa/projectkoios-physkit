# src/python/projectkoios/physkit/periodic/basis.py

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.lattice import (
    ReciprocalLattice1D,
    ReciprocalLattice3D,
)

type FloatArray = NDArray[np.float64]
type IntArray = NDArray[np.int64]
type ReciprocalIndex3D = tuple[int, int, int]

__all__ = ["ReciprocalModeBasis1D", "ReciprocalModeBasis3D"]


class ReciprocalModeBasis1D:
    """
    Reciprocal-lattice Fourier basis for a periodic one-dimensional field.

    The basis functions are

    .. math::

        \\phi_m(x)=e^{iG_mx},

    with

    .. math::

        G_m=mG_0.

    This basis is independent of the physical field being represented. It
    may be used for electronic wavefunctions, lattice-displacement fields,
    electromagnetic fields, or other periodic quantities.

    Parameters
    ----------
    reciprocal_lattice
        Reciprocal lattice defining ``G_0``.
    modes
        One-dimensional sequence of integer reciprocal-lattice modes.

    Attributes
    ----------
    reciprocal_lattice
        Reciprocal lattice associated with the basis.
    modes
        Read-only array of integer mode indices.
    vectors
        Read-only array of reciprocal vectors ``G_m``.
    """

    def __init__(
        self,
        reciprocal_lattice: ReciprocalLattice1D,
        modes: Sequence[int] | IntArray,
    ) -> None:
        self.reciprocal_lattice = reciprocal_lattice

        # Normalize the public input to a consistent internal representation.
        self.modes: IntArray = np.array(
            modes,
            dtype=np.int64,
            copy=True,
        )

        self.check_args()

        # Convert integer mode labels through the canonical reciprocal-lattice
        # operation so this basis shares the package's public lattice identity.
        self.vectors: FloatArray = self.reciprocal_lattice.vector(self.modes)

        # A basis must not change after dependent operators are constructed.
        self.modes.setflags(write=False)
        self.vectors.setflags(write=False)

    @classmethod
    def from_mode_limit(
        cls,
        reciprocal_lattice: ReciprocalLattice1D,
        mode_limit: int,
    ) -> ReciprocalModeBasis1D:
        """
        Construct a symmetric reciprocal-mode basis.

        Parameters
        ----------
        reciprocal_lattice
            Reciprocal lattice defining ``G_0``.
        mode_limit
            Maximum absolute integer mode ``M``.

        Returns
        -------
        ReciprocalModeBasis1D
            Basis containing modes ``-M`` through ``+M``.

        Raises
        ------
        TypeError
            If ``mode_limit`` is not an integer.
        ValueError
            If ``mode_limit`` is negative.
        """
        if not isinstance(
            mode_limit,
            (int, np.integer),
        ):
            raise TypeError("mode_limit must be an integer.")

        if mode_limit < 0:
            raise ValueError("mode_limit must be nonnegative.")

        modes: IntArray = np.arange(
            -mode_limit,
            mode_limit + 1,
            dtype=np.int64,
        )

        return cls(
            reciprocal_lattice=reciprocal_lattice,
            modes=modes,
        )

    @property
    def size(self) -> int:
        """
        Return the number of reciprocal modes.

        Returns
        -------
        int
            Number of basis functions.
        """
        return self.modes.size

    def check_args(self) -> None:
        """
        Validate the reciprocal-mode basis.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the reciprocal lattice is not the canonical one-dimensional identity.
        ValueError
            If the modes are empty, multidimensional, or duplicated.
        """
        if type(self.reciprocal_lattice) is not ReciprocalLattice1D:
            raise TypeError("reciprocal_lattice must be ReciprocalLattice1D")

        if self.modes.ndim != 1:
            raise ValueError("modes must be one-dimensional.")

        if self.modes.size == 0:
            raise ValueError("modes cannot be empty.")

        if np.unique(self.modes).size != self.modes.size:
            raise ValueError("modes must be unique.")

    def shifted_wavevectors(
        self,
        wavevector: float,
    ) -> FloatArray:
        """
        Calculate shifted reciprocal vectors.

        Parameters
        ----------
        wavevector
            Bloch wavevector ``k`` or phonon wavevector ``q``.

        Returns
        -------
        FloatArray
            Shifted vectors ``wavevector + G_m``.

        Raises
        ------
        TypeError
            If ``wavevector`` is not scalar.
        """
        if not np.isscalar(wavevector):
            raise TypeError("wavevector must be a scalar.")

        return self.vectors + float(wavevector)


class ReciprocalModeBasis3D:
    """Represent an ordered finite basis of 3D reciprocal-lattice modes.

    Parameters
    ----------
    reciprocal_lattice
        Canonical three-dimensional reciprocal lattice.
    modes
        Nonempty ordered inventory of unique integer triplets.

    Attributes
    ----------
    reciprocal_lattice
        Exact reciprocal-lattice owner retained by the basis.
    modes
        Immutable integer array with shape ``(number_of_modes, 3)``.
    vectors
        Immutable Cartesian reciprocal vectors with the same shape.

    Notes
    -----
    Explicit construction preserves caller order. The symmetric-cutoff factory uses
    first-index outer, second-index middle, and third-index inner ordering.
    """

    def __init__(
        self,
        reciprocal_lattice: ReciprocalLattice3D,
        modes: Sequence[ReciprocalIndex3D] | IntArray,
    ) -> None:
        raw_modes = np.asarray(modes)
        if raw_modes.size == 0:
            raise ValueError("modes cannot be empty")
        if raw_modes.ndim != 2 or raw_modes.shape[1] != 3:
            raise ValueError("modes must have shape (number_of_modes, 3)")
        if raw_modes.dtype.kind not in "iu":
            raise TypeError("modes must contain integers excluding booleans")
        if raw_modes.dtype.kind == "u" and np.any(raw_modes > np.iinfo(np.int64).max):
            raise ValueError("modes must fit signed 64-bit integer storage")

        self.reciprocal_lattice = reciprocal_lattice
        self.modes: IntArray = np.array(raw_modes, dtype=np.int64, copy=True)
        self.check_args()
        self.vectors: FloatArray = self.reciprocal_lattice.vector(self.modes)
        self.modes.setflags(write=False)
        self.vectors.setflags(write=False)

    @classmethod
    def from_mode_limit(
        cls,
        reciprocal_lattice: ReciprocalLattice3D,
        mode_limit: int,
    ) -> ReciprocalModeBasis3D:
        """Construct a symmetric cubic reciprocal-mode inventory.

        Parameters
        ----------
        reciprocal_lattice
            Canonical three-dimensional reciprocal lattice.
        mode_limit
            Maximum absolute integer component in every direction.

        Returns
        -------
        ReciprocalModeBasis3D
            Basis containing every triplet in ``[-M, M]^3``.

        Raises
        ------
        TypeError
            If ``mode_limit`` is not a built-in integer.
        ValueError
            If ``mode_limit`` is negative.
        """
        if type(mode_limit) is not int:
            raise TypeError("mode_limit must be a built-in integer")
        if mode_limit < 0:
            raise ValueError("mode_limit must be nonnegative")
        modes = np.asarray(
            [
                (first, second, third)
                for first in range(-mode_limit, mode_limit + 1)
                for second in range(-mode_limit, mode_limit + 1)
                for third in range(-mode_limit, mode_limit + 1)
            ],
            dtype=np.int64,
        )
        return cls(reciprocal_lattice=reciprocal_lattice, modes=modes)

    @property
    def size(self) -> int:
        """Return the represented basis dimension."""
        return int(self.modes.shape[0])

    def check_args(self) -> None:
        """Validate nominal lattice ownership and unique mode identity.

        Raises
        ------
        TypeError
            If the reciprocal lattice is not the canonical 3D identity.
        ValueError
            If duplicate reciprocal indices are present.
        """
        if type(self.reciprocal_lattice) is not ReciprocalLattice3D:
            raise TypeError("reciprocal_lattice must be ReciprocalLattice3D")
        if np.unique(self.modes, axis=0).shape[0] != self.modes.shape[0]:
            raise ValueError("modes must be unique")

    def shifted_wavevectors(
        self,
        wavevector: Sequence[float] | FloatArray,
    ) -> FloatArray:
        """Add one finite Cartesian wave vector to every reciprocal vector.

        Parameters
        ----------
        wavevector
            Finite Cartesian wave vector with shape ``(3,)``.

        Returns
        -------
        FloatArray
            Shifted vectors with shape ``(number_of_modes, 3)``.

        Raises
        ------
        TypeError
            If ``wavevector`` does not contain real numeric values.
        ValueError
            If ``wavevector`` does not have shape ``(3,)`` or is nonfinite.
        """
        raw_wavevector = np.asarray(wavevector)
        if raw_wavevector.dtype.kind not in "iuf":
            raise TypeError("wavevector must contain real numeric values")
        normalized = np.asarray(raw_wavevector, dtype=np.float64)
        if normalized.shape != (3,):
            raise ValueError("wavevector must have shape (3,)")
        if not np.all(np.isfinite(normalized)):
            raise ValueError("wavevector must contain finite values")
        return self.vectors + normalized
