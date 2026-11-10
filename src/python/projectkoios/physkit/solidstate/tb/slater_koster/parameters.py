"""Explicit shared parameter groups for homopolar ``sp3s*`` bonds."""

from __future__ import annotations

from dataclasses import dataclass

from .sp3s_star import HomopolarSp3sStarTwoCenterIntegrals


@dataclass(frozen=True, slots=True, kw_only=True)
class HomopolarSp3sStarParameterGroup:
    """Identify one explicitly shared seven-integral parameter set."""

    identifier: str
    integrals: HomopolarSp3sStarTwoCenterIntegrals

    def __post_init__(self) -> None:
        if type(self.identifier) is not str:
            raise TypeError("identifier must be a built-in string")
        if not self.identifier:
            raise ValueError("identifier must be nonempty")
        if type(self.integrals) is not HomopolarSp3sStarTwoCenterIntegrals:
            raise TypeError("integrals must be HomopolarSp3sStarTwoCenterIntegrals")
