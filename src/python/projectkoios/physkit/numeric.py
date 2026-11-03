# physkit/numeric.py
# Eugene Joseph M. Ragasa, 2026

from __future__ import annotations

import numpy as np

from projectkoios.physkit.types import ArrayLike


def as_f64_array(x: ArrayLike) -> np.ndarray:
    """
    Coerce input to a float64 NumPy array.
    """
    return np.asarray(x, dtype=np.float64)
