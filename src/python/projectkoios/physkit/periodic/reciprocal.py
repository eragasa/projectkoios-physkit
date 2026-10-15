# src/python/projectkoios/physkit/periodic/reciprocal.py

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.lattice import DirectLattice1D


class ReciprocalLattice1D:
    """
    One-dimensional reciprocal Bravais lattice.

    The direct-lattice primitive basis is

    .. math::

        \\mathcal A
        =
        \\{a_1\\},

    and the reciprocal-lattice primitive basis is

    .. math::

        \\mathcal B
        =
        \\{b_1\\},

    where

    .. math::

        a_1b_1
        =
        2\\pi.

    Therefore,

    .. math::

        b_1
        =
        \\frac{2\\pi}{a_1}.

    Every reciprocal-lattice vector is

    .. math::

        G_m
        =
        mb_1,

    where ``m`` is an integer.

    Parameters
    ----------
    direct_lattice:
        Direct lattice from which the reciprocal lattice is constructed.

    Attributes
    ----------
    direct_lattice:
        Associated one-dimensional direct lattice.
    primitive_basis:
        Primitive reciprocal-lattice basis vector ``b_1``. In one dimension,
        the vector is represented as a scalar.
    """

    def __init__(
        self,
        direct_lattice: DirectLattice1D,
    ) -> None:
        self.direct_lattice = direct_lattice

        self.check_args()

        # The reciprocal primitive basis is derived from the direct primitive
        # basis. It must not be supplied independently because the two bases
        # are constrained by a_1 b_1 = 2 pi.
        self.primitive_basis: float = 2.0 * np.pi / self.direct_lattice.lattice_constant

    def check_args(self) -> None:
        """
        Validate the reciprocal-lattice arguments.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If ``direct_lattice`` is not a ``DirectLattice1D`` instance.
        """
        if not isinstance(
            self.direct_lattice,
            DirectLattice1D,
        ):
            raise TypeError("direct_lattice must be a DirectLattice1D instance.")

    def vector(
        self,
        mode: int,
    ) -> float:
        """
        Construct a reciprocal-lattice vector.

        Parameters
        ----------
        mode:
            Integer reciprocal-lattice mode ``m``.

        Returns
        -------
        float
            Reciprocal-lattice vector ``G_m = mb_1``.

        Raises
        ------
        TypeError
            If ``mode`` is not an integer.
        """
        if not isinstance(
            mode,
            (int, np.integer),
        ):
            raise TypeError("mode must be an integer.")

        return float(mode) * self.primitive_basis
