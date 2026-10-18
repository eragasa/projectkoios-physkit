"""Public façade checks for unit-aware finite-periodic matrix operations."""

from projectkoios.physkit.solidstate import tb
from projectkoios.physkit.solidstate.tb import finite_periodic


def test_tb_facade_exports_exact_finite_periodic_identities() -> None:
    names = (
        "FinitePeriodicMatrixFourierEvaluation",
        "FinitePeriodicMatrixFourierSeries",
        "FinitePeriodicMatrixFourierTensorEvaluator",
        "FinitePeriodicMatrixFourierTransformer",
        "FinitePeriodicMatrixFourierTransformResult",
        "FinitePeriodicMatrixHermiticityAnalyzer",
        "FinitePeriodicMatrixHermiticityDiagnostics",
        "FinitePeriodicMatrixParsevalAnalyzer",
        "FinitePeriodicMatrixParsevalDiagnostics",
        "FinitePeriodicMatrixSamples",
    )

    for name in names:
        assert getattr(tb, name) is getattr(finite_periodic, name)
