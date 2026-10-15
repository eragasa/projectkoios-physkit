# src/python/projectkoios/physkit/qm/models/base.py

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TypeVar


class BaseQuantumModel(ABC):  # noqa: B024 - nominal abstract ownership root
    """Base physical specification of a quantum-mechanical system.

    A model defines the physical system. It does not define a numerical
    representation, solution method, or computational result.
    """


class QuantumResult(ABC):  # noqa: B024 - nominal abstract result root
    """
    Base result produced by a quantum-mechanical solver.
    """

    # Specific result families inherit from this class.
    pass


# This placeholder represents the exact result type returned by a
# particular solver.
QuantumResultType = TypeVar(
    "QuantumResultType",
    bound=QuantumResult,
    covariant=True,
)


class QuantumSolver[QuantumResultType: QuantumResult](
    ABC,
):
    """
    Base interface for a quantum-mechanical solver.
    """

    @abstractmethod
    def solve(self) -> QuantumResultType:
        """
        Solve the physical model and return the result.
        """

        # Each concrete solver supplies its own mathematics.
        raise NotImplementedError
