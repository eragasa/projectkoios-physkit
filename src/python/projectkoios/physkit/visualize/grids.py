from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from projectkoios.physkit.core.grids import (
    CartesianGrid1D,
    CartesianGrid2D,
    CartesianGrid3D,
)


def plot_grid_1d(grid: CartesianGrid1D) -> tuple[Figure, Axes]:
    """
    Plot the points of a 1D Cartesian grid.
    """
    fig, ax = plt.subplots()

    y = np.zeros_like(grid.x)

    ax.scatter(grid.x, y)
    ax.set_xlabel("x")
    ax.set_yticks([])
    ax.set_title("CartesianGrid1D")
    ax.grid(True)

    return fig, ax


def plot_grid_2d(
    grid: CartesianGrid2D, *, show_lines: bool = True
) -> tuple[Figure, Axes]:
    """
    Plot the points of a 2D Cartesian grid.
    """
    fig, ax = plt.subplots()

    X, Y = grid.mesh

    if show_lines:
        for i in range(grid.Nx):
            ax.plot(X[i, :], Y[i, :], linewidth=0.8)

        for j in range(grid.Ny):
            ax.plot(X[:, j], Y[:, j], linewidth=0.8)

    ax.scatter(X.ravel(), Y.ravel(), s=20)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("CartesianGrid2D")
    ax.set_aspect("equal", adjustable="box")

    return fig, ax


def plot_grid_3d(grid: CartesianGrid3D, *, max_points: int = 5000):
    """
    Plot the points of a 3D Cartesian grid.

    For large grids, this function downsamples the point cloud.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    points = grid.points

    if points.shape[0] > max_points:
        idx = np.linspace(0, points.shape[0] - 1, max_points).astype(int)
        points = points[idx]

    ax.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        s=5,
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("CartesianGrid3D")

    return fig, ax
