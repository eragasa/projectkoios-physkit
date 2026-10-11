# src/python/projectkoios/physkit/solidstate/electronic/models.py

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.lattice import DirectLattice1D
from projectkoios.physkit.potentials.periodic1d.cosine import (
    CosinePotential1D,
)
from projectkoios.physkit.units import ScalarQuantity, Unitless

type ComplexArray = NDArray[np.complex128]


class PeriodicPotentialFourier1D(ABC):
    """
    Abstract Fourier representation of a periodic potential.
    """

    @abstractmethod
    def coefficients(
        self,
        transferred_wavevectors,
    ) -> ComplexArray:
        """
        Evaluate potential Fourier coefficients.

        Parameters
        ----------
        transferred_wavevectors:
            Reciprocal transfers ``G' - G``.

        Returns
        -------
        ComplexArray
            Fourier coefficients ``V_(G' - G)``.
        """
        raise NotImplementedError


class CosinePotentialFourier1D(
    PeriodicPotentialFourier1D
):
    """
    Fourier representation of a cosine potential.

    The real-space potential is

    .. math::

        V(x)=V_0\\cos(G_0x).

    Its only nonzero Fourier coefficients are

    .. math::

        V_{+G_0}=V_{-G_0}=\\frac{V_0}{2}.

    Parameters
    ----------
    amplitude:
        Cosine amplitude ``V_0``.
    lattice:
        Direct lattice defining ``G_0=2*pi/a``.
    """

    def __init__(
        self,
        amplitude: float,
        lattice: DirectLattice1D,
    ) -> None:
        self.amplitude: float = float(amplitude)
        self.lattice = lattice

    def coefficients(
        self,
        transferred_wavevectors,
    ) -> ComplexArray:
        """
        Evaluate the cosine-potential Fourier coefficients.

        Parameters
        ----------
        transferred_wavevectors:
            Reciprocal transfers ``G' - G``.

        Returns
        -------
        ComplexArray
            Fourier coefficients with the same shape as the input.
        """
        transfers = np.asarray(
            transferred_wavevectors,
            dtype=np.float64,
        )

        fundamental_vector = (
            2.0
            * np.pi
            / self.lattice.lattice_constant
        )

        flattened_transfers = transfers.ravel()
        if flattened_transfers.size == 0:
            return np.zeros(transfers.shape, dtype=np.complex128)
        nearest_modes = np.rint(
            flattened_transfers / fundamental_vector
        ).astype(np.int64)
        reciprocal_modes = nearest_modes * fundamental_vector
        representable = np.isclose(flattened_transfers, reciprocal_modes)
        mode_indices = tuple((int(mode),) for mode in nearest_modes)
        periodic_potential = CosinePotential1D(
            amplitude=ScalarQuantity(self.amplitude, Unitless())
        )
        coefficients = periodic_potential.fourier_coefficients(
            mode_indices=mode_indices
        ).coefficients.magnitude.copy()
        coefficients[~representable] = 0.0
        return coefficients.reshape(transfers.shape)


class PlaneWaveBlochModel1D:
    """
    Physical model for a periodic one-electron problem.

    Parameters
    ----------
    lattice:
        One-dimensional direct lattice.
    potential:
        Fourier representation of the periodic potential.
    mass:
        Electron or effective particle mass.
    hbar:
        Reduced Planck constant in the selected unit system.
    """

    def __init__(
        self,
        lattice: DirectLattice1D,
        potential: PeriodicPotentialFourier1D,
        *,
        mass: float = 1.0,
        hbar: float = 1.0,
    ) -> None:
        self.lattice = lattice
        self.potential = potential
        self.mass: float = float(mass)
        self.hbar: float = float(hbar)

        self.check_args()

    def check_args(self) -> None:
        """
        Validate the electronic model.

        Returns
        -------
        None
        """
        if self.mass <= 0.0:
            raise ValueError(
                "mass must be positive."
            )

        if self.hbar <= 0.0:
            raise ValueError(
                "hbar must be positive."
            )
