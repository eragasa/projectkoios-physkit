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
    """Return a conjugate-supported sparse chain with analytic Bloch matrices."""
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


def test__SparseBlockHoppingBlochEvaluator__result__retains_sparse_matrices() -> None:
    """Retain sparse result, model, and matrix ownership through evaluation."""
    model = _model()

    result = SparseBlockHoppingBlochEvaluator().evaluate(
        model,
        MatrixQuantity(np.asarray([[0.0], [0.25]]), Unitless()),
    )

    # Compare the correlated ownership contract as one value rather than as an
    # assertion train over the result and each returned matrix.
    assert (
        type(result),
        result.model,
        tuple(type(matrix) for matrix in result.bloch_matrices),
    ) == (
        SparseBlockHoppingBlochEvaluation,
        model,
        (ComplexSparseMatrixQuantity, ComplexSparseMatrixQuantity),
    )


def test__SparseBlockHoppingBlochEvaluator__gamma__matches_analytic_matrix() -> None:
    """Match the analytic sparse Bloch matrix at the zone center."""
    result = SparseBlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.asarray([[0.0]]), Unitless()),
    )

    # At q=0 both directed hoppings contribute with phase one.
    np.testing.assert_array_equal(
        result.bloch_matrices[0].to_csr().toarray(),
        np.asarray([[2.0, 1.0], [1.0, 3.0]], dtype=np.complex128),
    )


def test__SparseBlockHoppingBlochEvaluator__quarter_turn__matches_matrix() -> None:
    """Match the analytic sparse Bloch matrix at one quarter reciprocal turn."""
    result = SparseBlockHoppingBlochEvaluator().evaluate(
        _model(),
        MatrixQuantity(np.asarray([[0.25]]), Unitless()),
    )

    # The ±1 displacements acquire ±i phases, distinguishing the exponent sign
    # while retaining the Hermitian conjugate pairing.
    np.testing.assert_allclose(
        result.bloch_matrices[0].to_csr().toarray(),
        np.asarray([[2.0, 1.0j], [-1.0j, 3.0]], dtype=np.complex128),
        rtol=0.0,
        atol=1.0e-15,
    )


def test__SparseBlockHoppingBlochEvaluator__empty_grid__retains_correlation() -> None:
    """An empty reciprocal grid returns an empty exact sparse result inventory."""
    wave_vectors = MatrixQuantity(np.empty((0, 1), dtype=np.float64), Unitless())

    result = SparseBlockHoppingBlochEvaluator().evaluate(_model(), wave_vectors)

    # No placeholder zero matrix is created when there are no reciprocal points.
    assert result.bloch_matrices == ()
