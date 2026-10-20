"""Unit-aware physical adapters for finite-periodic matrix Fourier operations.

ProjectKoios Numerics owns centered meshes, Fourier arithmetic, tensor-product
contractions, Parseval reductions, and cyclic Hermiticity pairing. This package owns
the physical correlation that those generic kernels deliberately omit: one explicit
energy unit, energy reference, and ordered-basis identity across samples,
coefficients, evaluations, and diagnostics.

The records and actions choose no support, cutoff, tolerance, or scientific acceptance
policy. Reduced coordinates and integer modes are dimensionless; matrix samples,
coefficients, evaluated matrices, and Hermiticity defects carry the declared energy
unit. Parseval quantities carry its square.
"""

from .actions import (
    FinitePeriodicMatrixFourierTensorEvaluator,
    FinitePeriodicMatrixFourierTransformer,
)
from .diagnostics import (
    FinitePeriodicMatrixHermiticityAnalyzer,
    FinitePeriodicMatrixHermiticityDiagnostics,
    FinitePeriodicMatrixParsevalAnalyzer,
    FinitePeriodicMatrixParsevalDiagnostics,
)
from .records import (
    FinitePeriodicMatrixFourierEvaluation,
    FinitePeriodicMatrixFourierSeries,
    FinitePeriodicMatrixFourierTransformResult,
    FinitePeriodicMatrixSamples,
)

__all__ = [
    "FinitePeriodicMatrixFourierEvaluation",
    "FinitePeriodicMatrixFourierSeries",
    "FinitePeriodicMatrixFourierTensorEvaluator",
    "FinitePeriodicMatrixFourierTransformResult",
    "FinitePeriodicMatrixFourierTransformer",
    "FinitePeriodicMatrixHermiticityAnalyzer",
    "FinitePeriodicMatrixHermiticityDiagnostics",
    "FinitePeriodicMatrixParsevalAnalyzer",
    "FinitePeriodicMatrixParsevalDiagnostics",
    "FinitePeriodicMatrixSamples",
]
