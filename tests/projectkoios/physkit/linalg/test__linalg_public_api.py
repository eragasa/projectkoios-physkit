"""Public façade tests for generic linear-algebra identities."""

import projectkoios.physkit.linalg as linalg
import projectkoios.physkit.linalg.matrices as matrices
import projectkoios.physkit.linalg.matrix.norm as matrix_norm


def test__linalg_facade__retains_defining_module_identities() -> None:
    """The package façade exports exact objects rather than wrappers."""
    assert linalg.CsrSparseMatrix is matrices.CsrSparseMatrix
    assert linalg.MatrixInput is matrices.MatrixInput
    assert linalg.MatrixVerifier is matrices.MatrixVerifier
    assert linalg.VerifiedMatrix is matrices.VerifiedMatrix
    assert linalg.VerifiedMatrixUtils is matrices.VerifiedMatrixUtils
    assert linalg.MatrixNorm is matrix_norm.MatrixNorm
    assert linalg.MatrixNormEnum is matrix_norm.MatrixNormEnum
    assert (
        linalg.ScipyCsrSpectralNormControls is matrix_norm.ScipyCsrSpectralNormControls
    )
    assert (
        linalg.SparseSpectralNormOptionsEnum
        is matrix_norm.SparseSpectralNormOptionsEnum
    )
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
