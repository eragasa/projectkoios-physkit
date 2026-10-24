"""Public façade for eager represented spectral decompositions."""

from __future__ import annotations

from projectkoios.physkit.linalg.eigen.base import SpectralDecompositionABC
from projectkoios.physkit.linalg.eigen.eigensystem import Eigensystem
from projectkoios.physkit.linalg.eigen.factory import SpectralDecompositionFactory
from projectkoios.physkit.linalg.eigen.types import SpectralMatrixQuantity


class SpectralDecomposition:
    """Eagerly compute and expose one supported represented decomposition."""

    __slots__ = ("_implementation",)

    def __init__(self, matrix: SpectralMatrixQuantity) -> None:
        """Construct the selected product and complete its maximum-rank solve."""
        self._implementation: SpectralDecompositionABC = (
            SpectralDecompositionFactory().create(matrix)
        )

    @property
    def matrix(self) -> SpectralMatrixQuantity:
        """Return the exact matrix retained by the selected implementation."""
        return self._implementation.matrix

    @property
    def max_rank(self) -> int:
        """Return matrix dimension N, not a numerical matrix-rank estimate."""
        return self._implementation.max_rank

    def eigensystem(
        self,
        n_rank: int | None = None,
    ) -> Eigensystem:
        """Return the retained full eigensystem or its lowest prefix."""
        return self._implementation.eigensystem(n_rank)
