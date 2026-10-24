"""Matrix-quantity families admitted by represented spectral contracts."""

from projectkoios.physkit.units import (
    ComplexMatrixQuantity,
    ComplexSparseMatrixQuantity,
    MatrixQuantity,
    SparseMatrixQuantity,
)

type RealSpectralMatrixQuantity = MatrixQuantity | SparseMatrixQuantity
type ComplexSpectralMatrixQuantity = ComplexMatrixQuantity | ComplexSparseMatrixQuantity
type SpectralMatrixQuantity = RealSpectralMatrixQuantity | ComplexSpectralMatrixQuantity
