"""Public API tests for matrix operation packages."""

import projectkoios.physkit.linalg.matrix as matrix
import projectkoios.physkit.linalg.matrix.norm as matrix_norm


def test__matrix_public_api__exposes_exact_norm_subpackage() -> None:
    """The matrix package exposes its norm family without flattening it."""
    assert matrix.norm is matrix_norm
    assert matrix.__all__ == ["norm"]
