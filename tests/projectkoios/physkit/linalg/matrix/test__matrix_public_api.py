"""Public API tests for the structural matrix-operations package."""

import projectkoios.physkit.linalg.matrix as matrix


def test__matrix_public_api__does_not_flatten_operation_identities() -> None:
    """The structural package exports no subpackage-owned objects."""
    assert matrix.__all__ == []
    assert not hasattr(matrix, "MatrixNorm")
    assert not hasattr(matrix, "MatrixNormEnum")
    assert not hasattr(matrix, "SparseSpectralNormOptionsEnum")
