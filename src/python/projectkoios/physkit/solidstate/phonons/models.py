# src/python/projectkoios/physkit/solidstate/phonons/models.py

from __future__ import annotations

import numpy as np

from projectkoios.physkit.periodic.lattice import DirectLattice1D


class MonatomicChainModel1D:
    """
    One-dimensional monatomic harmonic chain.

    Parameters
    ----------
    lattice:
        Direct lattice with nearest-neighbor spacing ``a``.
    atomic_mass:
        Mass ``M`` of each atom.
    spring_constant:
        Nearest-neighbor spring constant ``K``.
    """

    def __init__(
        self,
        lattice: DirectLattice1D,
        *,
        atomic_mass: float,
        spring_constant: float,
    ) -> None:
        self.lattice = lattice
        self.atomic_mass: float = float(atomic_mass)
        self.spring_constant: float = float(spring_constant)

        self.check_args()

    def check_args(self) -> None:
        """
        Validate the monatomic-chain model.

        Returns
        -------
        None
        """
        if not np.isfinite(self.atomic_mass):
            raise ValueError("atomic_mass must be finite.")

        if self.atomic_mass <= 0.0:
            raise ValueError("atomic_mass must be positive.")

        if not np.isfinite(self.spring_constant):
            raise ValueError("spring_constant must be finite.")

        if self.spring_constant <= 0.0:
            raise ValueError("spring_constant must be positive.")
