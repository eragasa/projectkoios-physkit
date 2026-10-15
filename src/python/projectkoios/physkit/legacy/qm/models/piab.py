"""Deprecated particle-in-a-box model interfaces.

Use :mod:`projectkoios.physkit.qm.piab1d` for the maintained one-dimensional model.
"""

from __future__ import annotations

from abc import ABC

from projectkoios.physkit.qm.models.base import (
    BaseQuantumModel,
    QuantumResult,
    QuantumSolver,
)


class Piab(
    BaseQuantumModel,
    ABC,
):
    """
    Base physical specification of a particle-in-a-box model.

    This class identifies the physical model family. Dimensional
    geometry and particle parameters are defined by subclasses.
    """

    # Piab1D, Piab2D, and Piab3D inherit from this class.
    pass


class PiabResult(
    QuantumResult,
    ABC,
):
    """
    Base result produced by a particle-in-a-box solver.
    """

    # Dimension-specific PIAB result classes inherit from this class.
    pass


class PiabSolver[PiabResultType: PiabResult](
    QuantumSolver[PiabResultType],
    ABC,
):
    """
    Base interface for a particle-in-a-box solver.
    """

    # The solve() method is inherited from QuantumSolver.
    pass
