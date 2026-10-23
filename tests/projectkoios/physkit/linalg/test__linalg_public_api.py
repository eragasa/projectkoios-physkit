"""Public façade tests for generic linear-algebra identities."""

import projectkoios.physkit.linalg as linalg
import projectkoios.physkit.linalg.eigenproblems as eigenproblems
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
    eigenproblem_names = (
        "RealSymmetricEigendecomposer",
        "RealSymmetricEigendecomposition",
        "RealSymmetricOperator",
    )

    for name in matrix_names:
        assert getattr(linalg, name) is getattr(matrices, name)
    for name in eigenproblem_names:
        assert getattr(linalg, name) is getattr(eigenproblems, name)
    assert linalg.__all__ == sorted((*matrix_names, *eigenproblem_names))
