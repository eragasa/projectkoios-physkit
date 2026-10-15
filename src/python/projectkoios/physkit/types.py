# physkit/types.py

from __future__ import annotations

import numpy as np

# Scalar or numeric NumPy array (most common scientific pattern)
Scalar = float
FloatArray = np.ndarray
ArrayLike = Scalar | FloatArray
