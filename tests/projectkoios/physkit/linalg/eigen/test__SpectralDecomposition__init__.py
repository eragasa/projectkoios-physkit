"""Verification of the public eager spectral-decomposition façade."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.linalg.eigen import SpectralDecomposition
from projectkoios.physkit.units import MatrixQuantity, SparseMatrixQuantity, Unitless


def test__init__eagerly_retains_complete_eigensystem_and_matrix_identity() -> None:
    matrix = MatrixQuantity(np.diag([1.0, 2.0]), Unitless())

    decomposition = SpectralDecomposition(matrix)

    assert decomposition.matrix is matrix
    assert decomposition.max_rank == 2
    assert decomposition.eigensystem().matrix is matrix
    assert decomposition.eigensystem().is_complete
    assert not hasattr(decomposition, "decompose")


def test__init__retains_complete_sparse_tridiagonal_boundary() -> None:
    matrix = SparseMatrixQuantity.from_csr(
        sparse.csr_array(
            np.array(
                [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]],
                dtype=np.float64,
            )
        ),
        Unitless(),
    )

    with pytest.raises(ValueError, match="tridiagonal"):
        SpectralDecomposition(matrix)
