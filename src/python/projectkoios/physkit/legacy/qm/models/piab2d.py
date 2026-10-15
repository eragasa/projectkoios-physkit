"""Deprecated, incomplete two-dimensional PIAB declarations."""

from __future__ import annotations

from abc import ABC

from .piab import Piab, PiabResult, PiabSolver


class Piab2D(Piab):
    """
    Physical specification of a two-dimensional particle in a box.
    """

    pass


class Piab2DResult(PiabResult):
    """
    Base result for a two-dimensional particle-in-a-box model.
    """

    pass


class Piab2DSolver[Piab2DResultT: Piab2DResult](
    PiabSolver[Piab2DResultT],
    ABC,
):
    """
    Base solver for a two-dimensional particle-in-a-box model.
    """

    pass


class Piab2DAnalyticalResult(Piab2DResult):
    pass


class Piab2DAnalyticalSolver(
    Piab2DSolver[Piab2DAnalyticalResult],
):
    @property
    def result(self) -> Piab2DAnalyticalResult:
        raise NotImplementedError


class Piab2DNumericalResult(Piab2DResult):
    pass


class Piab2DNumericalSolver(
    Piab2DSolver[Piab2DNumericalResult],
):
    @property
    def result(self) -> Piab2DNumericalResult:
        raise NotImplementedError
