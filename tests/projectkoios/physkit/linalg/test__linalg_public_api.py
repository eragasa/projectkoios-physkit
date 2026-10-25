"""Public façade tests for generic linear-algebra identities."""

import projectkoios.physkit.linalg as linalg
import projectkoios.physkit.linalg.matrices as matrices
import projectkoios.physkit.linalg.matrix_norms as matrix_norms


def test__linalg_facade__retains_defining_module_identities() -> None:
    """The package façade exports exact objects rather than wrappers."""
    matrix_names = (
        "CsrSparseMatrix",
        "MatrixInput",
        "MatrixVerifier",
        "VerifiedMatrix",
        "VerifiedMatrixUtils",
    )
    matrix_norm_names = (
        "MatrixNorm",
        "MatrixNormEnum",
        "ScipyCsrSpectralNormControls",
        "SparseSpectralNormOptionsEnum",
    )
    for name in matrix_names:
        assert getattr(linalg, name) is getattr(matrices, name)
    for name in matrix_norm_names:
        assert getattr(linalg, name) is getattr(matrix_norms, name)
    assert linalg.__all__ == [
        "CsrSparseMatrix",
        "MatrixInput",
        "MatrixNorm",
        "MatrixNormEnum",
        "MatrixVerifier",
        "ScipyCsrSpectralNormControls",
        "SparseSpectralNormOptionsEnum",
        "VerifiedMatrix",
        "VerifiedMatrixUtils",
    ]
