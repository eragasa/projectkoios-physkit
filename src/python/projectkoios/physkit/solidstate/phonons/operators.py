# src/python/projectkoios/physkit/solidstate/phonons/operators.py

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.solidstate.phonons.models import (
    MonatomicChainModel1D,
)

type ComplexArray = NDArray[np.complex128]


class MonatomicDynamicalMatrix1D:
    """
    Dynamical matrix for a monatomic nearest-neighbor chain.

    For a monatomic one-dimensional chain, the dynamical matrix is a
    one-by-one matrix:

    .. math::

        D(q)
        =
        \\frac{2K}{M}
        \\left[1-\\cos(qa)\\right].

    Its eigenvalue is the squared phonon angular frequency:

    .. math::

        D(q)e(q)
        =
        \\omega^2(q)e(q).
    """

    def __init__(
        self,
        model: MonatomicChainModel1D,
    ) -> None:
        self.model = model

    def construct_matrix(
        self,
        phonon_wavevector: float,
    ) -> ComplexArray:
        """
        Construct the dynamical matrix at one phonon wavevector.

        Parameters
        ----------
        phonon_wavevector:
            Phonon wavevector ``q``.

        Returns
        -------
        ComplexArray
            One-by-one Hermitian dynamical matrix.
        """
        q = float(phonon_wavevector)
        a = self.model.lattice.lattice_constant
        mass = self.model.atomic_mass
        spring_constant = self.model.spring_constant

        eigenvalue = 2.0 * spring_constant / mass * (1.0 - np.cos(q * a))

        return np.array(
            [[eigenvalue]],
            dtype=np.complex128,
        )
