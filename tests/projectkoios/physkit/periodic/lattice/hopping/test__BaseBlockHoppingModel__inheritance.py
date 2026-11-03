"""Nominal inheritance tests for dense and sparse block hopping models."""

import pytest

from projectkoios.physkit.periodic.lattice.hopping import (
    BaseBlockHoppingModel,
    BaseBlockHoppingTerm,
    DenseBlockHoppingModel,
    DenseBlockHoppingTerm,
    SparseBlockHoppingModel,
    SparseBlockHoppingTerm,
)


def test__inheritance__dense_and_sparse_terms__share_nominal_base() -> None:
    """Both coefficient representations have one nominal term owner."""
    assert DenseBlockHoppingTerm.__bases__ == (BaseBlockHoppingTerm,)
    assert SparseBlockHoppingTerm.__bases__ == (BaseBlockHoppingTerm,)


def test__inheritance__dense_and_sparse_models__share_nominal_base() -> None:
    """Both model representations have one nominal model owner."""
    assert DenseBlockHoppingModel.__bases__ == (BaseBlockHoppingModel,)
    assert SparseBlockHoppingModel.__bases__ == (BaseBlockHoppingModel,)


def test__abstract_bases__own_shared_behavior_and_reject_construction() -> None:
    """Abstract owners retain shared derivation and cannot exist incomplete."""
    assert DenseBlockHoppingTerm.basis_size is BaseBlockHoppingTerm.basis_size
    assert SparseBlockHoppingTerm.basis_size is BaseBlockHoppingTerm.basis_size
    assert DenseBlockHoppingModel.basis_size is BaseBlockHoppingModel.basis_size
    assert SparseBlockHoppingModel.basis_size is BaseBlockHoppingModel.basis_size
    with pytest.raises(TypeError, match="abstract"):
        BaseBlockHoppingTerm()  # type: ignore[abstract]
    with pytest.raises(TypeError, match="abstract"):
        BaseBlockHoppingModel()  # type: ignore[abstract]
