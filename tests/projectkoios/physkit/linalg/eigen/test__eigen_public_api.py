"""Public façade tests for represented spectral decomposition."""

import projectkoios.physkit.linalg.eigen as eigen
import projectkoios.physkit.linalg.eigen.base as base
import projectkoios.physkit.linalg.eigen.eigensystem as eigensystem_module
import projectkoios.physkit.linalg.eigen.errors as errors
import projectkoios.physkit.linalg.eigen.facade as facade
import projectkoios.physkit.linalg.eigen.types as eigen_types


def test__eigen_facade__retains_supported_object_identities() -> None:
    """The package exports exact contracts without concrete implementation names."""
    type_names = (
        "ComplexSpectralMatrixQuantity",
        "RealSpectralMatrixQuantity",
        "SpectralMatrixQuantity",
    )

    for name in type_names:
        assert getattr(eigen, name) is getattr(eigen_types, name)
    assert eigen.SpectralDecompositionABC is base.SpectralDecompositionABC
    assert eigen.Eigensystem is eigensystem_module.Eigensystem
    assert eigen.SpectralDecomposition is facade.SpectralDecomposition
    assert (
        eigen.UnsupportedSpectralDecompositionError
        is errors.UnsupportedSpectralDecompositionError
    )
    assert not hasattr(eigen, "SpectralDecompositionFactory")
    assert not hasattr(eigen, "RealSymmetricSpectralDecomposition")
    assert eigen.__all__ == [
        "ComplexSpectralMatrixQuantity",
        "Eigensystem",
        "RealSpectralMatrixQuantity",
        "SpectralDecomposition",
        "SpectralDecompositionABC",
        "SpectralMatrixQuantity",
        "UnsupportedSpectralDecompositionError",
    ]
