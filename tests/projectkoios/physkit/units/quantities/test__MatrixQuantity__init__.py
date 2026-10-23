"""Verification of immutable ``MatrixQuantity`` storage."""

import numpy as np
import pytest

from projectkoios.physkit.units.quantities import MatrixQuantity, PhysicalUnit


def test_owns_immutable_binary64_values() -> None:
    source = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
    quantity = MatrixQuantity(source, PhysicalUnit("joule"))
    source[:] = 0.0

    np.testing.assert_array_equal(quantity.magnitude, [[1.0, 2.0], [3.0, 4.0]])
    assert quantity.magnitude.dtype == np.dtype("<f8")
    assert not quantity.magnitude.flags.writeable


@pytest.mark.parametrize("dtype", [np.int64, np.float32])
def test__MatrixQuantity__init__accepts_real_numeric_arrays(
    dtype: type[np.generic],
) -> None:
    source = np.array([[1, 2], [3, 4]], dtype=dtype)

    quantity = MatrixQuantity(source, PhysicalUnit("joule"))

    np.testing.assert_array_equal(
        quantity.magnitude,
        [[1.0, 2.0], [3.0, 4.0]],
    )
    assert quantity.magnitude.dtype == np.dtype("<f8")


@pytest.mark.parametrize(
    "source",
    [
        np.array([[True, False]], dtype=np.bool_),
        np.array([["1.0", "2.0"]]),
        np.array([[b"1.0", b"2.0"]]),
        np.array([[1.0, 2.0]], dtype=object),
        np.array([[1.0 + 0.0j, 2.0 + 1.0j]], dtype=np.complex128),
    ],
)
def test__MatrixQuantity__init__rejects_nonreal_numeric_arrays(
    source: np.ndarray,
) -> None:
    with pytest.raises(TypeError, match="real numeric values excluding booleans"):
        MatrixQuantity(source, PhysicalUnit("joule"))
