"""Construction tests for sparse matrix-valued hopping terms."""

import numpy as np
import pytest
from scipy import sparse

from projectkoios.physkit.periodic.lattice.finite_domain import (
    LatticeDimension,
    LatticeDisplacement,
)
from projectkoios.physkit.periodic.lattice.hopping import SparseBlockHoppingTerm
from projectkoios.physkit.units import ComplexSparseMatrixQuantity, Unitless


def test__init__canonical_sparse_square_block__retains_sparse_ownership() -> None:
    """A canonical sparse square coefficient retains its CSR representation."""
    coefficient = ComplexSparseMatrixQuantity.from_csr(
        sparse.csr_array(np.diag([1.0 + 2.0j, 3.0 - 4.0j])),
        Unitless(),
    )
    term = SparseBlockHoppingTerm(
        LatticeDisplacement(LatticeDimension.ONE, (1,)),
        coefficient,
    )

    assert term.coefficient is coefficient
    assert term.basis_size == 2
    assert term.energy_unit == Unitless()
    assert term.storage_kind == "sparse"
    assert term.coefficient.nonzero_count == 2


def test__init__nonsquare_or_empty_shape__rejects() -> None:
    """Sparse hopping blocks require a nonempty square represented space."""
    for shape in ((1, 2), (0, 0)):
        coefficient = ComplexSparseMatrixQuantity.from_csr(
            sparse.csr_array(shape, dtype=np.complex128),
            Unitless(),
        )
        with pytest.raises(ValueError, match="nonempty square"):
            SparseBlockHoppingTerm(
                LatticeDisplacement(LatticeDimension.ONE, (0,)),
                coefficient,
            )
