# src/python/projectkoios/physkit/periodic/basis.py

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.reciprocal import ReciprocalLattice1D

type FloatArray = NDArray[np.float64]
type IntArray = NDArray[np.int64]


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
    reciprocal_lattice:
        Reciprocal lattice defining ``G_0``.
    modes:
        One-dimensional sequence of integer reciprocal-lattice modes.

    Attributes
    ----------
    reciprocal_lattice:
        Reciprocal lattice associated with the basis.
    modes:
        Read-only array of integer mode indices.
    vectors:
        Read-only array of reciprocal vectors ``G_m``.
    """

    def __init__(
        self,
        reciprocal_lattice: ReciprocalLattice1D,
        modes: Sequence[int] | IntArray,
    ) -> None:
        self.reciprocal_lattice = reciprocal_lattice

        # Normalize the public input to a consistent internal representation.
        self.modes: IntArray = np.asarray(
            modes,
            dtype=np.int64,
        )

        self.check_args()

        # Convert integer mode labels to physical reciprocal vectors.
        self.vectors: FloatArray = (
            self.modes.astype(
                np.float64,
                copy=False,
            )
            * self.reciprocal_lattice.fundamental_vector
        )

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
        reciprocal_lattice:
            Reciprocal lattice defining ``G_0``.
        mode_limit:
            Maximum absolute integer mode ``M``.

        Returns
        -------
        ReciprocalModeBasis1D
            Basis containing modes ``-M`` through ``+M``.
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
        ValueError
            If the modes are empty, multidimensional, or duplicated.
        """
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
        wavevector:
            Bloch wavevector ``k`` or phonon wavevector ``q``.

        Returns
        -------
        FloatArray
            Shifted vectors ``wavevector + G_m``.
        """
        if not np.isscalar(wavevector):
            raise TypeError("wavevector must be a scalar.")

        return self.vectors + float(wavevector)
