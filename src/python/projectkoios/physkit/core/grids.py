# src/python/projectkoios/physkit/grids.py
# 2026, Eugene Joseph M. Ragasa

from __future__ import annotations

from numbers import Real

import numpy as np
from numpy.typing import NDArray

type FloatArray = NDArray[np.float64]


class CartesianAxis:
    """
    Uniformly sampled one-dimensional Cartesian coordinate axis.

    The axis represents the interval from ``lower`` to ``upper`` using
    ``N`` uniformly spaced coordinate samples.

    Parameters
    ----------
    name : str
        Name or label of the coordinate axis, such as ``"x"``, ``"y"``,
        or ``"z"``.
    lower : float
        Lower bound of the coordinate interval.
    upper : float
        Upper bound of the coordinate interval. It must be greater than
        ``lower``.
    N : int
        Number of coordinate samples. It must be at least 2.
    endpoint : bool, optional
        Whether ``upper`` is included in the coordinate samples.
        The default is ``True``.

    Attributes
    ----------
    name : str
        Name or label of the coordinate axis.
    lower : float
        Lower bound of the coordinate interval.
    upper : float
        Upper bound of the coordinate interval.
    N : int
        Number of coordinate samples.
    endpoint : bool
        Whether the upper bound is included in the coordinate samples.
    L : float
        Length of the coordinate interval, equal to ``upper - lower``.
    d : float
        Uniform spacing between adjacent coordinate samples.
    values : numpy.ndarray
        One-dimensional array containing the coordinate samples. The
        array has shape ``(N,)`` and dtype ``numpy.float64``.

    Raises
    ------
    TypeError
        If ``name`` is not a string, ``lower`` or ``upper`` is not a
        real scalar, ``N`` is not an integer, or ``endpoint`` is not a
        Boolean.
    ValueError
        If either bound is nonfinite, ``upper`` is not greater than
        ``lower``, or ``N`` is less than 2.

    Notes
    -----
    The interval length is

    .. math::

        L = x_{\\mathrm{upper}} - x_{\\mathrm{lower}}.

    When ``endpoint=True``, the upper bound is included and the sample
    spacing is

    .. math::

        d = \\frac{L}{N-1}.

    When ``endpoint=False``, the upper bound is excluded and the sample
    spacing is

    .. math::

        d = \\frac{L}{N}.

    The endpoint-excluded convention is suitable for sampling periodic
    domains without duplicating the first sample at the upper boundary.

    Examples
    --------
    Construct a closed axis containing both bounds:

    >>> axis = CartesianAxis("x", 0.0, 1.0, 5)
    >>> axis.values
    array([0.  , 0.25, 0.5 , 0.75, 1.  ])
    >>> axis.d
    0.25

    Construct an endpoint-excluded axis:

    >>> axis = CartesianAxis("x", 0.0, 1.0, 4, endpoint=False)
    >>> axis.values
    array([0.  , 0.25, 0.5 , 0.75])
    """

    def __init__(
        self,
        name: str,
        lower: float,
        upper: float,
        N: int,
        *,
        endpoint: bool = True,
    ) -> None:
        if not isinstance(name, str):
            raise TypeError("name must be a string.")

        if isinstance(lower, bool) or not isinstance(lower, Real):
            raise TypeError("lower must be a real scalar.")

        if isinstance(upper, bool) or not isinstance(upper, Real):
            raise TypeError("upper must be a real scalar.")

        if isinstance(N, bool) or not isinstance(N, (int, np.integer)):
            raise TypeError(f"N{name} must be an integer.")

        if not isinstance(endpoint, bool):
            raise TypeError(f"endpoint_{name} must be a bool.")

        self.name: str = name
        self.lower: float = float(lower)
        self.upper: float = float(upper)
        self.N: int = int(N)
        self.endpoint: bool = endpoint

        self.check_args()

        self.L: float = self.upper - self.lower
        self.length: float = self.L

        self.values: FloatArray = np.linspace(
            start=self.lower,
            stop=self.upper,
            num=self.N,
            endpoint=self.endpoint,
            dtype=np.float64,
        )
        self.values.setflags(write=False)

        self.d: float = float(self.values[1] - self.values[0])
        self.delta: float = self.d

    def check_args(self) -> None:
        """
        Validate the numerical axis parameters.

        Raises
        ------
        ValueError
            If either coordinate bound is nonfinite, ``upper`` is not
            greater than ``lower``, or ``N`` is less than 2.
        """
        if not np.isfinite(self.lower):
            raise ValueError(f"{self.name}_min must be finite.")

        if not np.isfinite(self.upper):
            raise ValueError(f"{self.name}_max must be finite.")

        if self.upper <= self.lower:
            raise ValueError(f"{self.name}_max must be greater than {self.name}_min.")

        if self.N < 2:
            raise ValueError(f"N{self.name} must be at least 2.")

    def __repr__(self) -> str:
        """
        Return an unambiguous string representation of the axis.
        """
        return (
            f"CartesianAxis("
            f"name={self.name!r}, "
            f"lower={self.lower}, "
            f"upper={self.upper}, "
            f"N={self.N}, "
            f"endpoint={self.endpoint})"
        )


class CartesianGrid1D:
    """
    One-dimensional uniform Cartesian grid.

    The grid defines the discrete coordinate representation used for
    fields that depend on the Cartesian coordinate ``x``.

    Parameters
    ----------
    x_min : float
        Lower boundary of the spatial domain.
    x_max : float
        Upper boundary of the spatial domain.
    Nx : int
        Number of coordinate samples along the x-axis.
    endpoint : bool, optional
        Whether ``x_max`` is included as a coordinate sample. The
        default is ``True``.

    Attributes
    ----------
    x : numpy.ndarray
        Read-only Cartesian x-coordinate samples.
    x_axis : CartesianAxis
        Cartesian x-axis metadata.
    shape : tuple[int]
        Shape of a scalar field represented on the grid.
    size : int
        Total number of grid points.

    Examples
    --------
    >>> grid = CartesianGrid1D(
    ...     x_min=0.0,
    ...     x_max=1.0,
    ...     Nx=5,
    ... )
    >>> grid.x
    array([0.  , 0.25, 0.5 , 0.75, 1.  ])
    >>> grid.dx
    0.25
    >>> grid.shape
    (5,)
    """

    def __init__(
        self,
        x_min: float,
        x_max: float,
        Nx: int,
        *,
        endpoint: bool = True,
    ) -> None:
        self.x_axis = CartesianAxis(
            name="x",
            lower=x_min,
            upper=x_max,
            N=Nx,
            endpoint=endpoint,
        )

        self.x_min: float = self.x_axis.lower
        self.x_max: float = self.x_axis.upper
        self.Nx: int = self.x_axis.N
        self.endpoint: bool = self.x_axis.endpoint
        self.x: FloatArray = self.x_axis.values
        self.Lx: float = self.x_axis.L
        self.dx: float = self.x_axis.d
        self.shape: tuple[int] = (self.Nx,)
        self.size: int = self.Nx
        self.delta: float = self.dx

    def __repr__(self) -> str:
        return (
            f"CartesianGrid1D("
            f"x_min={self.x_min}, "
            f"x_max={self.x_max}, "
            f"Nx={self.Nx}, "
            f"endpoint={self.endpoint})"
        )


class ActiveSet1D:
    """
    Active grid indices on a one-dimensional Cartesian grid.

    Parameters
    ----------
    grid : CartesianGrid1D
        Grid on which the active indices are defined.
    indices : numpy.ndarray
        One-dimensional array of unique active grid indices.

    Attributes
    ----------
    grid : CartesianGrid1D
        Grid associated with the active indices.
    indices : numpy.ndarray
        Active grid indices with dtype ``numpy.int64``.

    Raises
    ------
    ValueError
        If ``indices`` is not one-dimensional, contains repeated
        indices, or contains an index outside the grid.
    """

    def __init__(
        self,
        grid: CartesianGrid1D,
        indices: NDArray[np.int64],
    ) -> None:
        self.grid: CartesianGrid1D = grid

        self.indices: NDArray[np.int64] = np.asarray(
            indices,
            dtype=np.int64,
        )

        if self.indices.ndim != 1:
            raise ValueError("indices must be one-dimensional.")

        if np.any(self.indices < 0):
            raise ValueError("indices must not contain negative values.")

        if np.any(self.indices >= self.grid.Nx):
            raise ValueError("indices must be less than grid.Nx.")

        if np.unique(self.indices).size != self.indices.size:
            raise ValueError("indices must not contain repeated values.")

    @property
    def coordinates(self) -> FloatArray:
        """
        Return the coordinates at the active grid indices.

        Returns
        -------
        numpy.ndarray
            Active coordinate values.
        """
        return self.grid.x[self.indices]

    @property
    def size(self) -> int:
        """
        Return the number of active grid indices.

        Returns
        -------
        int
            Number of active indices.
        """
        return int(self.indices.size)

    def __repr__(self) -> str:
        return f"ActiveSet1D(size={self.size}, indices={self.indices!r})"


class CartesianGrid2D:
    """
    Two-dimensional uniform Cartesian grid.

    The grid represents scalar fields sampled over Cartesian coordinates
    ``x`` and ``y``. A scalar field on this grid has shape ``(Nx, Ny)``.

    Parameters
    ----------
    x_min : float
        Lower boundary of the x-coordinate interval.
    x_max : float
        Upper boundary of the x-coordinate interval.
    Nx : int
        Number of samples along the x-axis.
    y_min : float
        Lower boundary of the y-coordinate interval.
    y_max : float
        Upper boundary of the y-coordinate interval.
    Ny : int
        Number of samples along the y-axis.
    endpoint_x : bool, optional
        Whether ``x_max`` is included. The default is ``True``.
    endpoint_y : bool, optional
        Whether ``y_max`` is included. The default is ``True``.

    Attributes
    ----------
    x : numpy.ndarray
        Read-only Cartesian x-coordinate samples.
    y : numpy.ndarray
        Read-only Cartesian y-coordinate samples.
    shape : tuple[int, int]
        Shape of a scalar field represented on the grid.
    size : int
        Total number of grid points.
    delta : tuple[float, float]
        Grid spacings ``(dx, dy)``.
    """

    def __init__(
        self,
        x_min: float,
        x_max: float,
        Nx: int,
        y_min: float,
        y_max: float,
        Ny: int,
        *,
        endpoint_x: bool = True,
        endpoint_y: bool = True,
    ) -> None:
        self.x_axis = CartesianAxis(
            name="x",
            lower=x_min,
            upper=x_max,
            N=Nx,
            endpoint=endpoint_x,
        )
        self.y_axis = CartesianAxis(
            name="y",
            lower=y_min,
            upper=y_max,
            N=Ny,
            endpoint=endpoint_y,
        )

        self.x_min: float = self.x_axis.lower
        self.x_max: float = self.x_axis.upper
        self.Nx: int = self.x_axis.N
        self.endpoint_x: bool = self.x_axis.endpoint
        self.y_min: float = self.y_axis.lower
        self.y_max: float = self.y_axis.upper
        self.Ny: int = self.y_axis.N
        self.endpoint_y: bool = self.y_axis.endpoint
        self.x: FloatArray = self.x_axis.values
        self.y: FloatArray = self.y_axis.values
        self.Lx: float = self.x_axis.L
        self.Ly: float = self.y_axis.L
        self.dx: float = self.x_axis.d
        self.dy: float = self.y_axis.d
        self.shape: tuple[int, int] = (self.Nx, self.Ny)
        self.size: int = self.Nx * self.Ny
        self.delta: tuple[float, float] = (self.dx, self.dy)

    @property
    def mesh(
        self,
    ) -> tuple[FloatArray, FloatArray]:
        """
        Return the Cartesian coordinate mesh.

        Returns
        -------
        tuple[FloatArray, FloatArray]
            Coordinate arrays ``X`` and ``Y``, each with shape
            ``(Nx, Ny)``.
        """
        return np.meshgrid(
            self.x,
            self.y,
            indexing="ij",
        )

    @property
    def points(self) -> FloatArray:
        """
        Return all Cartesian grid points.

        Returns
        -------
        FloatArray
            Coordinate pairs with shape ``(Nx * Ny, 2)``.
        """
        X, Y = self.mesh

        return np.column_stack(
            (
                X.ravel(),
                Y.ravel(),
            )
        )

    def __repr__(self) -> str:
        """Return an unambiguous representation of the grid."""
        return (
            f"CartesianGrid2D("
            f"x_min={self.x_min}, "
            f"x_max={self.x_max}, "
            f"Nx={self.Nx}, "
            f"y_min={self.y_min}, "
            f"y_max={self.y_max}, "
            f"Ny={self.Ny}, "
            f"endpoint_x={self.endpoint_x}, "
            f"endpoint_y={self.endpoint_y})"
        )


class CartesianGrid3D:
    """
    Three-dimensional uniform Cartesian grid.

    The grid represents scalar fields sampled over Cartesian coordinates
    ``x``, ``y``, and ``z``. A scalar field on this grid has shape
    ``(Nx, Ny, Nz)``.

    Parameters
    ----------
    x_min : float
        Lower boundary of the x-coordinate interval.
    x_max : float
        Upper boundary of the x-coordinate interval.
    Nx : int
        Number of samples along the x-axis.
    y_min : float
        Lower boundary of the y-coordinate interval.
    y_max : float
        Upper boundary of the y-coordinate interval.
    Ny : int
        Number of samples along the y-axis.
    z_min : float
        Lower boundary of the z-coordinate interval.
    z_max : float
        Upper boundary of the z-coordinate interval.
    Nz : int
        Number of samples along the z-axis.
    endpoint_x : bool, optional
        Whether ``x_max`` is included. The default is ``True``.
    endpoint_y : bool, optional
        Whether ``y_max`` is included. The default is ``True``.
    endpoint_z : bool, optional
        Whether ``z_max`` is included. The default is ``True``.

    Attributes
    ----------
    x : numpy.ndarray
        Read-only Cartesian x-coordinate samples.
    y : numpy.ndarray
        Read-only Cartesian y-coordinate samples.
    z : numpy.ndarray
        Read-only Cartesian z-coordinate samples.
    shape : tuple[int, int, int]
        Shape of a scalar field represented on the grid.
    size : int
        Total number of grid points.
    delta : tuple[float, float, float]
        Grid spacings ``(dx, dy, dz)``.
    """

    def __init__(
        self,
        x_min: float,
        x_max: float,
        Nx: int,
        y_min: float,
        y_max: float,
        Ny: int,
        z_min: float,
        z_max: float,
        Nz: int,
        *,
        endpoint_x: bool = True,
        endpoint_y: bool = True,
        endpoint_z: bool = True,
    ) -> None:
        self.x_axis = CartesianAxis(
            name="x",
            lower=x_min,
            upper=x_max,
            N=Nx,
            endpoint=endpoint_x,
        )
        self.y_axis = CartesianAxis(
            name="y",
            lower=y_min,
            upper=y_max,
            N=Ny,
            endpoint=endpoint_y,
        )
        self.z_axis = CartesianAxis(
            name="z",
            lower=z_min,
            upper=z_max,
            N=Nz,
            endpoint=endpoint_z,
        )

        self.x_min: float = self.x_axis.lower
        self.x_max: float = self.x_axis.upper
        self.Nx: int = self.x_axis.N
        self.endpoint_x: bool = self.x_axis.endpoint
        self.y_min: float = self.y_axis.lower
        self.y_max: float = self.y_axis.upper
        self.Ny: int = self.y_axis.N
        self.endpoint_y: bool = self.y_axis.endpoint
        self.z_min: float = self.z_axis.lower
        self.z_max: float = self.z_axis.upper
        self.Nz: int = self.z_axis.N
        self.endpoint_z: bool = self.z_axis.endpoint
        self.x: FloatArray = self.x_axis.values
        self.y: FloatArray = self.y_axis.values
        self.z: FloatArray = self.z_axis.values
        self.Lx: float = self.x_axis.L
        self.Ly: float = self.y_axis.L
        self.Lz: float = self.z_axis.L
        self.dx: float = self.x_axis.d
        self.dy: float = self.y_axis.d
        self.dz: float = self.z_axis.d
        self.shape: tuple[int, int, int] = (self.Nx, self.Ny, self.Nz)
        self.size: int = self.Nx * self.Ny * self.Nz
        self.delta: tuple[float, float, float] = (
            self.dx,
            self.dy,
            self.dz,
        )

    @property
    def mesh(
        self,
    ) -> tuple[FloatArray, FloatArray, FloatArray]:
        """
        Return the Cartesian coordinate mesh.

        Returns
        -------
        tuple[FloatArray, FloatArray, FloatArray]
            Coordinate arrays ``X``, ``Y``, and ``Z``, each with shape
            ``(Nx, Ny, Nz)``.
        """
        return np.meshgrid(
            self.x,
            self.y,
            self.z,
            indexing="ij",
        )

    @property
    def points(self) -> FloatArray:
        """
        Return all Cartesian grid points.

        Returns
        -------
        FloatArray
            Coordinate triples with shape ``(Nx * Ny * Nz, 3)``.
        """
        X, Y, Z = self.mesh

        return np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

    def __repr__(self) -> str:
        """Return an unambiguous representation of the grid."""
        return (
            f"CartesianGrid3D("
            f"x_min={self.x_min}, "
            f"x_max={self.x_max}, "
            f"Nx={self.Nx}, "
            f"y_min={self.y_min}, "
            f"y_max={self.y_max}, "
            f"Ny={self.Ny}, "
            f"z_min={self.z_min}, "
            f"z_max={self.z_max}, "
            f"Nz={self.Nz}, "
            f"endpoint_x={self.endpoint_x}, "
            f"endpoint_y={self.endpoint_y}, "
            f"endpoint_z={self.endpoint_z})"
        )
