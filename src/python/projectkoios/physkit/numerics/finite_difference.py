# physkit/numerics/finitedifference.py
# Eugene Ragasa, 2026

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix

from projectkoios.physkit.core.boundaries import (
    BoundaryCondition,
)
from projectkoios.physkit.core.boundaries import (
    DirichletBoundaryCondition as DirichletBC,
)
from projectkoios.physkit.core.grids import CartesianGrid1D

_DEFAULT_DIRICHLET_BC = DirichletBC(value=0.0)


class Laplacian1D:
    """
    Finite-difference matrix representation of the 1D Cartesian Laplacian.

    In one Cartesian dimension,

        ∇² = d²/dx²

    This class builds matrix representations of d²/dx² on a CartesianGrid1D
    with specified boundary conditions.

    Important
    ---------
    Homogeneous Dirichlet boundary conditions remove the boundary point from
    the set of unknowns.

    For example, if the grid has Nx points and both endpoints are Dirichlet,
    the returned matrix has shape

        (Nx - 2, Nx - 2)

    because the boundary values are known.

    Homogeneous Neumann boundary conditions keep the boundary point as an
    unknown and modify the boundary stencil using ghost-point reflection.

    Periodic boundary conditions keep all points and connect the first and
    last grid points. Periodic grids should use endpoint=False.
    """

    def __init__(
        self,
        grid: CartesianGrid1D,
        *,
        bc_x_min: BoundaryCondition = _DEFAULT_DIRICHLET_BC,
        bc_x_max: BoundaryCondition = _DEFAULT_DIRICHLET_BC,
    ) -> None:
        self.grid: CartesianGrid1D = grid
        self.bc_x_min: BoundaryCondition = bc_x_min
        self.bc_x_max: BoundaryCondition = bc_x_max
        self.check_args()
        self.check_boundary_conditions()

        self.x_active_min_idx = self.get_x_active_min_idx()
        self.x_active_max_idx = self.get_x_active_max_idx()
        self.x_active = self.get_x_active()
        self.Nx_active = self.x_active.shape[0]
        self.shape: tuple[int, int] = (self.Nx_active, self.Nx_active)

        self.matrix: csr_matrix = self.build()

    def check_args(self) -> None:
        if not isinstance(self.grid, CartesianGrid1D):
            raise TypeError("grid must be a CartesianGrid1D.")
        if not isinstance(self.bc_x_min, BoundaryCondition):
            raise TypeError("bc_x_min must be a BoundaryCondition")
        if not isinstance(self.bc_x_max, BoundaryCondition):
            raise TypeError("bc_x_max must be a BoundaryCondition")

    def check_boundary_conditions(self) -> None:
        supported_bcs: tuple[type[BoundaryCondition], ...] = (DirichletBC,)

        if not isinstance(self.bc_x_min, supported_bcs):
            raise TypeError(f"bc_x_min:{self.bc_x_min} is not a supported BC")
        if not isinstance(self.bc_x_max, supported_bcs):
            raise TypeError(f"bc_x_min:{self.bc_x_max} is not a supported BC")

        if isinstance(self.bc_x_min, DirichletBC) and self.bc_x_min.value != 0:
            raise ValueError("Only DirichletBC(value=0.0) is supported at x_min.")

        if isinstance(self.bc_x_max, DirichletBC) and self.bc_x_max.value != 0:
            raise ValueError("Only DirichletBC(value=0.0) is supported at x_max.")

    def get_x_active_min_idx(self) -> int:
        x_active_min_idx = 0
        if isinstance(self.bc_x_min, DirichletBC):
            x_active_min_idx = 1
        return x_active_min_idx

    def get_x_active_max_idx(self) -> int:
        x_active_max_idx = self.grid.Nx
        if isinstance(self.bc_x_max, DirichletBC):
            x_active_max_idx = self.grid.Nx - 1

        return x_active_max_idx

    def get_x_active(self) -> np.ndarray:
        return self.grid.x[self.x_active_min_idx : self.x_active_max_idx]

    def build(self) -> csr_matrix:
        """
        Build the homogeneous Dirichlet 1D Laplacian matrix.

        The returned matrix acts on the active unknown vector:

            u_active = [u_1, u_2, ..., u_{N-2}]^T

        for boundary values

            u_0 = 0
            u_{N-1} = 0
        """

        main = -2.0 * np.ones(self.Nx_active)
        off = np.ones(self.Nx_active - 1)

        L = sp.diags(
            diagonals=[off, main, off],
            offsets=[-1, 0, 1],
            shape=self.shape,
            format="csr",
        )
        L = L / self.grid.dx**2
        return L
