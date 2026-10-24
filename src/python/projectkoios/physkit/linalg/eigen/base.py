"""Nominal contracts for eager represented spectral decompositions."""

from __future__ import annotations

from abc import ABC, abstractmethod

from projectkoios.physkit.linalg.eigen.eigensystem import Eigensystem
from projectkoios.physkit.linalg.eigen.types import SpectralMatrixQuantity


class SpectralDecompositionABC(ABC):
    """Require an eagerly computed maximum-rank represented decomposition."""

    @property
    @abstractmethod
    def matrix(self) -> SpectralMatrixQuantity:
        """Return the exact matrix retained by the decomposition."""

    @property
    @abstractmethod
    def max_rank(self) -> int:
        """Return matrix dimension N, not a numerical matrix-rank estimate."""

    @abstractmethod
    def eigensystem(
        self,
        n_rank: int | None = None,
    ) -> Eigensystem:
        """Return the retained full eigensystem or its lowest prefix."""
