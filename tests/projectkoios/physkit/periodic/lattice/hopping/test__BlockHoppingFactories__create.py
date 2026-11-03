"""Dense and sparse block-hopping factory tests."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    BlockHoppingModel,
    BlockHoppingTerm,
    DenseBlockHoppingModel,
    DenseBlockHoppingTerm,
    SparseBlockHoppingModel,
    SparseBlockHoppingTerm,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    PhysicalUnit,
)


class TestBlockHoppingFactories:
    """Verify exact construction without representation conversion."""

    unit = PhysicalUnit("electron_volt")
    displacement = LatticeDisplacement(LatticeDimension.ONE, (0,))

    def test__term_factory__exact_quantity__selects_exact_term(self) -> None:
        dense_quantity = ComplexMatrixQuantity(np.eye(2), self.unit)
        sparse_quantity = ComplexSparseMatrixQuantity.from_csr(
            sparse.eye(2, format="csr", dtype=np.complex128),
            self.unit,
        )

        factory = BlockHoppingTerm()
        dense = factory.create(self.displacement, dense_quantity)
        sparse_term = factory.create(self.displacement, sparse_quantity)

        assert type(dense) is DenseBlockHoppingTerm
        assert dense.coefficient is dense_quantity
        assert type(sparse_term) is SparseBlockHoppingTerm
        assert sparse_term.coefficient is sparse_quantity

    def test__model_factory__homogeneous_terms__selects_exact_model(self) -> None:
        dense_term = DenseBlockHoppingTerm(
            self.displacement,
            ComplexMatrixQuantity(np.eye(1), self.unit),
        )
        sparse_term = SparseBlockHoppingTerm(
            self.displacement,
            ComplexSparseMatrixQuantity.from_csr(
                sparse.eye(1, format="csr", dtype=np.complex128),
                self.unit,
            ),
        )
        arguments = {
            "identifier": "one-state-chain",
            "dimension": LatticeDimension.ONE,
            "energy_unit": self.unit,
            "energy_reference": "zero",
            "basis_identifier": "one-state-basis",
        }

        factory = BlockHoppingModel()
        dense = factory.create(terms=(dense_term,), **arguments)
        sparse_model = factory.create(terms=(sparse_term,), **arguments)

        assert type(dense) is DenseBlockHoppingModel
        assert type(sparse_model) is SparseBlockHoppingModel

    def test__factories__unsupported_or_mixed_representation__rejects(self) -> None:
        dense_term = DenseBlockHoppingTerm(
            self.displacement,
            ComplexMatrixQuantity(np.eye(1), self.unit),
        )
        sparse_term = SparseBlockHoppingTerm(
            LatticeDisplacement(LatticeDimension.ONE, (1,)),
            ComplexSparseMatrixQuantity.from_csr(
                sparse.eye(1, format="csr", dtype=np.complex128),
                self.unit,
            ),
        )
        term_factory = BlockHoppingTerm()
        model_factory = BlockHoppingModel()
        with pytest.raises(TypeError, match="coefficient must be"):
            term_factory.create(
                self.displacement,
                object(),  # type: ignore[arg-type]
            )
        with pytest.raises(ValueError, match="terms must be nonempty"):
            model_factory.create(
                "empty",
                LatticeDimension.ONE,
                (),
                self.unit,
                "zero",
                "basis",
            )
        with pytest.raises(TypeError, match="only DenseBlockHoppingTerm"):
            model_factory.create(
                "mixed",
                LatticeDimension.ONE,
                (dense_term, sparse_term),  # type: ignore[arg-type]
                self.unit,
                "zero",
                "basis",
            )
