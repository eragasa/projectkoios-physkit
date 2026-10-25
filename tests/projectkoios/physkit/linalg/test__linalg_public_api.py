"""Public API tests for the structural linear-algebra package."""

import projectkoios.physkit.linalg as linalg


def test__linalg_public_api__does_not_flatten_owned_identities() -> None:
    """The structural package exports no module or subpackage-owned objects."""
    assert linalg.__all__ == []
    assert not hasattr(linalg, "MatrixVerifier")
    assert not hasattr(linalg, "MatrixNorm")
    assert not hasattr(linalg, "MatrixNormEnum")
    assert not hasattr(linalg, "ScipyCsrSpectralNormControls")
    assert not hasattr(linalg, "SparseSpectralNormOptionsEnum")
