"""Direct sparse Bloch-evaluation tests."""

import numpy as np
from scipy import sparse

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    SparseBlockHoppingModel,
    SparseBlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bloch import (
    SparseBlockHoppingBlochEvaluation,
    SparseBlockHoppingBlochEvaluator,
)
from projectkoios.physkit.units import (
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)


def _model() -> SparseBlockHoppingModel:
    unit = PhysicalUnit("electron_volt")
    blocks = (
        (-1, sparse.csr_array([[0.0, 0.0], [1.0, 0.0]], dtype=np.complex128)),
        (0, sparse.csr_array([[2.0, 0.0], [0.0, 3.0]], dtype=np.complex128)),
        (1, sparse.csr_array([[0.0, 1.0], [0.0, 0.0]], dtype=np.complex128)),
    )
    return SparseBlockHoppingModel(
        identifier="sparse-chain",
        dimension=LatticeDimension.ONE,
        terms=tuple(
            SparseBlockHoppingTerm(
                LatticeDisplacement(LatticeDimension.ONE, (displacement,)),
                ComplexSparseMatrixQuantity.from_csr(block, unit),
            )
            for displacement, block in blocks
        ),
        energy_unit=unit,
        energy_reference="zero",
        basis_identifier="a,b",
    )


def test__evaluate__sparse_blocks__returns_sparse_bloch_matrices() -> None:
    """Sparse coefficients remain sparse through direct Bloch summation."""
    model = _model()
    wave_vectors = MatrixQuantity(np.asarray([[0.0], [0.25]]), Unitless())

    result = SparseBlockHoppingBlochEvaluator().evaluate(model, wave_vectors)

    assert type(result) is SparseBlockHoppingBlochEvaluation
    assert result.model is model
    assert all(
        type(matrix) is ComplexSparseMatrixQuantity
        for matrix in result.bloch_matrices
    )
    assert np.array_equal(
        result.bloch_matrices[0].to_csr().toarray(),
        np.asarray([[2.0, 1.0], [1.0, 3.0]], dtype=np.complex128),
    )
    assert np.allclose(
        result.bloch_matrices[1].to_csr().toarray(),
        np.asarray([[2.0, 1.0j], [-1.0j, 3.0]], dtype=np.complex128),
        rtol=0.0,
        atol=1.0e-15,
    )


def test__evaluate__empty_grid__retains_empty_sparse_correlation() -> None:
    """An empty reciprocal grid returns an empty exact sparse result inventory."""
    wave_vectors = MatrixQuantity(np.empty((0, 1), dtype=np.float64), Unitless())

    result = SparseBlockHoppingBlochEvaluator().evaluate(_model(), wave_vectors)

    assert result.bloch_matrices == ()
