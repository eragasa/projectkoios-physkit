"""Dense/sparse parity and basis covariance for block-hopping Bloch sums."""

from __future__ import annotations

from typing import Never

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import (
    DenseBlockHoppingModel,
    DenseBlockHoppingTerm,
    SparseBlockHoppingModel,
    SparseBlockHoppingTerm,
)
from projectkoios.physkit.solidstate.tb.bloch import (
    DenseBlockHoppingBlochEvaluator,
    SparseBlockHoppingBlochEvaluator,
)
from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    PhysicalUnit,
    Unitless,
)

pytestmark = pytest.mark.numerical_verification


class TestDenseSparseBlockHoppingBlochParity:
    """Exercise representation parity and retained-basis covariance."""

    unit = PhysicalUnit("electron_volt")
    coordinates = MatrixQuantity(
        np.asarray(
            [
                [-0.5, -0.5],
                [-0.25, 0.125],
                [0.0, 0.0],
                [0.2, -0.3],
                [0.5, 0.5],
            ],
            dtype=np.float64,
        ),
        Unitless(),
    )

    @classmethod
    def _blocks(cls) -> tuple[tuple[tuple[int, int], np.ndarray], ...]:
        """Return complex three-state blocks with conjugate displacement support."""
        onsite = np.asarray(
            [
                [1.2, 0.1 + 0.2j, 0.0],
                [0.1 - 0.2j, -0.4, 0.3],
                [0.0, 0.3, 2.1],
            ],
            dtype=np.complex128,
        )
        x_positive = np.asarray(
            [
                [0.0, 0.25j, 0.0],
                [0.1, 0.0, -0.2j],
                [0.0, 0.0, 0.15],
            ],
            dtype=np.complex128,
        )
        y_positive = np.asarray(
            [
                [0.05, 0.0, -0.12],
                [0.0, -0.08j, 0.0],
                [0.2j, 0.0, 0.0],
            ],
            dtype=np.complex128,
        )
        return (
            ((-1, 0), x_positive.conjugate().T),
            ((0, -1), y_positive.conjugate().T),
            ((0, 0), onsite),
            ((0, 1), y_positive),
            ((1, 0), x_positive),
        )

    @classmethod
    def _dense_model(
        cls,
        blocks: tuple[tuple[tuple[int, int], np.ndarray], ...] | None = None,
    ) -> DenseBlockHoppingModel:
        """Construct the dense representation of the selected hopping blocks."""
        selected_blocks = cls._blocks() if blocks is None else blocks
        return DenseBlockHoppingModel(
            identifier="dense-two-dimensional-three-state-control",
            dimension=LatticeDimension.TWO,
            terms=tuple(
                DenseBlockHoppingTerm(
                    LatticeDisplacement(LatticeDimension.TWO, displacement),
                    ComplexMatrixQuantity(block, cls.unit),
                )
                for displacement, block in selected_blocks
            ),
            energy_unit=cls.unit,
            energy_reference="control-zero",
            basis_identifier="ordered-three-state-cell-basis",
        )

    @classmethod
    def _sparse_model(cls) -> SparseBlockHoppingModel:
        """Construct a CSR representation of the unchanged hopping blocks."""
        return SparseBlockHoppingModel(
            identifier="sparse-two-dimensional-three-state-control",
            dimension=LatticeDimension.TWO,
            terms=tuple(
                SparseBlockHoppingTerm(
                    LatticeDisplacement(LatticeDimension.TWO, displacement),
                    ComplexSparseMatrixQuantity.from_csr(
                        sparse.csr_array(block),
                        cls.unit,
                    ),
                )
                for displacement, block in cls._blocks()
            ),
            energy_unit=cls.unit,
            energy_reference="control-zero",
            basis_identifier="ordered-three-state-cell-basis",
        )

    def test__SparseBlockHoppingBlochEvaluator__storage__does_not_densify(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Complete sparse evaluation without invoking coefficient conversion."""

        def reject_dense_conversion(
            _quantity: ComplexSparseMatrixQuantity,
        ) -> Never:
            raise AssertionError(
                "sparse Bloch evaluation must not densify coefficients"
            )

        # A failing conversion sentinel turns the no-densification requirement into
        # an executable contract rather than an implementation assumption.
        monkeypatch.setattr(
            ComplexSparseMatrixQuantity,
            "to_dense",
            reject_dense_conversion,
        )

        result = SparseBlockHoppingBlochEvaluator().evaluate(
            self._sparse_model(),
            self.coordinates,
        )

        assert len(result.bloch_matrices) == self.coordinates.magnitude.shape[0]

    def test__DenseSparseBlochEvaluators__representation_parity__agree(self) -> None:
        """Produce the same Bloch matrices from equivalent dense and CSR blocks."""
        dense = DenseBlockHoppingBlochEvaluator().evaluate(
            self._dense_model(),
            self.coordinates,
        )
        sparse_evaluation = SparseBlockHoppingBlochEvaluator().evaluate(
            self._sparse_model(),
            self.coordinates,
        )

        # Stack the full reciprocal inventory so parity is one numerical requirement,
        # not one assertion per point and matrix entry.
        dense_values = np.stack(
            tuple(matrix.magnitude for matrix in dense.bloch_matrices),
            axis=0,
        )
        sparse_values = np.stack(
            tuple(
                matrix.to_csr().toarray()
                for matrix in sparse_evaluation.bloch_matrices
            ),
            axis=0,
        )
        np.testing.assert_allclose(
            sparse_values,
            dense_values,
            rtol=0.0,
            atol=2.0e-15,
        )

    def test__SparseBlockHoppingBlochEvaluator__result__is_canonical(
        self,
    ) -> None:
        """Return correlated energy-unit matrices in canonical CSR form."""
        result = SparseBlockHoppingBlochEvaluator().evaluate(
            self._sparse_model(),
            self.coordinates,
        )

        # Coordinate, unit, and canonical-storage correlation is one result contract.
        assert (
            result.fractional_wave_vectors,
            tuple(matrix.unit for matrix in result.bloch_matrices),
            tuple(
                matrix.to_csr().has_canonical_format
                for matrix in result.bloch_matrices
            ),
        ) == (
            self.coordinates,
            (self.unit,) * self.coordinates.magnitude.shape[0],
            (True,) * self.coordinates.magnitude.shape[0],
        )

    def test__DenseBlockHoppingBlochEvaluator__covariance__commutes_with_unitary(
        self,
    ) -> None:
        """Covary every Bloch matrix under one retained-basis unitary."""
        # The normalized three-point Fourier matrix is deterministic, complex, and
        # nontrivial in every basis direction, making covariance failures visible.
        indices = np.arange(3, dtype=np.float64)
        unitary = np.exp(
            2.0j * np.pi * np.outer(indices, indices) / 3.0
        ) / np.sqrt(3.0)
        # Apply one fixed retained-space basis change to every real-space block;
        # this changes representation without changing the represented operator.
        transformed_blocks = tuple(
            (
                displacement,
                unitary.conjugate().T @ block @ unitary,
            )
            for displacement, block in self._blocks()
        )

        source = DenseBlockHoppingBlochEvaluator().evaluate(
            self._dense_model(),
            self.coordinates,
        )
        transformed = DenseBlockHoppingBlochEvaluator().evaluate(
            self._dense_model(transformed_blocks),
            self.coordinates,
        )

        # Bloch summation is linear, so compare the complete reciprocal inventory
        # against one independently transformed array.
        source_values = np.stack(
            tuple(matrix.magnitude for matrix in source.bloch_matrices),
            axis=0,
        )
        transformed_values = np.stack(
            tuple(matrix.magnitude for matrix in transformed.bloch_matrices),
            axis=0,
        )
        expected = np.einsum(
            "ab,qbc,cd->qad",
            unitary.conjugate().T,
            source_values,
            unitary,
        )
        np.testing.assert_allclose(
            transformed_values,
            expected,
            rtol=0.0,
            atol=5.0e-15,
        )
