from dataclasses import dataclass
from enum import Enum, auto

import numpy as np


class GridType1D(Enum):
    LEFT_CLOSED = auto()  # [a, b)
    RIGHT_CLOSED = auto()  # (a, b]
    OPEN = auto()  # (a, b) (abstract open interval)
    INTERIOR = auto()  # (a, b) with N interior points (Dirichlet-style)
    MIDPOINT = auto()  # cell-centered
    CLOSED = auto()  # [a, b]


@dataclass(frozen=True)
class Grid1D:
    a: float
    b: float
    N: int
    grid_type: GridType1D = GridType1D.CLOSED

    def __post_init__(self):
        if self.b <= self.a:
            raise ValueError("Require b > a")
        if self.N <= 0:
            raise ValueError("Require N > 0")
        if self.grid_type is GridType1D.CLOSED and self.N < 2:
            raise ValueError("For CLOSED grid, require N >=2")

    @property
    def L(self):
        return self.b - self.a

    @property
    def dx(self) -> float:
        if self.grid_type is GridType1D.CLOSED:
            return self.L / (self.N - 1)
        elif self.grid_type is GridType1D.OPEN or self.grid_type is GridType1D.INTERIOR:
            return self.L / (self.N + 1)
        else:
            return self.L / self.N

    @property
    def x_arr(self) -> np.ndarray:
        dx = self.dx

        if self.grid_type is GridType1D.LEFT_CLOSED:
            return self.a + dx * np.arange(self.N)

        elif (
            self.grid_type is GridType1D.RIGHT_CLOSED
            or self.grid_type is GridType1D.OPEN
            or self.grid_type is GridType1D.INTERIOR
        ):
            return self.a + dx * (np.arange(self.N) + 1)

        elif self.grid_type is GridType1D.MIDPOINT:
            return self.a + dx * (np.arange(self.N) + 0.5)

        elif self.grid_type is GridType1D.CLOSED:
            return np.linspace(self.a, self.b, self.N)

        else:
            raise RuntimeError(f"Unhandled {self.grid_type}")

    @property
    def X(self) -> np.ndarray:
        r"""
        Discrete configuration `X` (implemented as an array of points).
        $x \in X$
        """
        return self.x_arr
