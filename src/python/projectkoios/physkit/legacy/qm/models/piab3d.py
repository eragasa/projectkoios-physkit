"""Deprecated, incomplete three-dimensional PIAB declarations."""

from __future__ import annotations

from abc import ABC

from .piab import Piab, PiabResult, PiabSolver


class Piab3D(Piab):
    """
    Physical specification of a three-dimensional particle in a box.
    """

    pass


class Piab3DResult(PiabResult):
    """
    Base result for a three-dimensional particle-in-a-box model.
    """

    pass


class Piab3DSolver[Piab3DResultT: Piab3DResult](
    PiabSolver[Piab3DResultT],
    ABC,
):
    """
    Base solver for a three-dimensional particle-in-a-box model.
    """

    pass


class Piab3DAnalyticalResult(Piab3DResult):
    pass


class Piab3DAnalyticalSolver(
    Piab3DSolver[Piab3DAnalyticalResult],
):
    @property
    def result(self) -> Piab3DAnalyticalResult:
        raise NotImplementedError


class Piab3DNumericalResult(Piab3DResult):
    pass


class Piab3DNumericalSolver(
    Piab3DSolver[Piab3DNumericalResult],
):
    @property
    def result(self) -> Piab3DNumericalResult:
        raise NotImplementedError
