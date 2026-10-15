# src/python/projectkoios/physkit/units/protocols.py
# Author: Eugene Joseph M. Ragasa
from __future__ import annotations

from enum import IntEnum
from typing import Protocol, runtime_checkable

from projectkoios.physkit.types import ArrayLike


@runtime_checkable
class UnitQuantityProtocol(Protocol):
    """
    Structural interface for a physkit physical quantity.

    A conforming quantity must:
    - define a nested Units(IntEnum)
    - provide convert(from_, to)
    - be stateless
    """

    class Units(IntEnum): ...

    @staticmethod
    def convert(
        value: ArrayLike, units_from: IntEnum, units_to: IntEnum
    ) -> ArrayLike: ...
