"""Public represented spectral-decomposition contracts and façade."""

from .base import SpectralDecompositionABC
from .eigensystem import Eigensystem
from .errors import UnsupportedSpectralDecompositionError
from .facade import SpectralDecomposition
from .types import (
    ComplexSpectralMatrixQuantity,
    RealSpectralMatrixQuantity,
    SpectralMatrixQuantity,
)

__all__ = [
    "ComplexSpectralMatrixQuantity",
    "Eigensystem",
    "RealSpectralMatrixQuantity",
    "SpectralDecomposition",
    "SpectralDecompositionABC",
    "SpectralMatrixQuantity",
    "UnsupportedSpectralDecompositionError",
]
