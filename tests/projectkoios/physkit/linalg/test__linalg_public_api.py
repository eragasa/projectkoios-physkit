"""Public façade tests for generic linear-algebra identities."""

import projectkoios.physkit.linalg as linalg
import projectkoios.physkit.linalg.matrices as matrices


def test__linalg_facade__retains_defining_module_identities() -> None:
    """The package façade exports exact objects rather than wrappers."""
    matrix_names = (
        "CsrSparseMatrix",
        "MatrixInput",
        "MatrixVerifier",
        "VerifiedMatrix",
        "VerifiedMatrixUtils",
    )
    for name in matrix_names:
        assert getattr(linalg, name) is getattr(matrices, name)
    assert linalg.__all__ == list(matrix_names)
