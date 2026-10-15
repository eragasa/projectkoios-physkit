"""Wavevector paths through one-, two-, and three-dimensional zones."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from projectkoios.physkit.periodic.lattice import FirstBrillouinZone1D

type FloatArray = NDArray[np.float64]
type IntArray = NDArray[np.int64]


class KPointPath1D:
    """Ordered path through a one-dimensional first Brillouin zone."""

    def __init__(self, values: FloatArray) -> None:
        self.values: FloatArray = np.asarray(values, dtype=np.float64)
        self.check_args()
        self.values.setflags(write=False)

    @classmethod
    def from_first_brillouin_zone(
        cls,
        zone: FirstBrillouinZone1D,
        number_of_points: int,
        *,
        endpoint: bool = True,
    ) -> KPointPath1D:
        """Uniformly sample a one-dimensional first Brillouin zone."""
        if not isinstance(zone, FirstBrillouinZone1D):
            raise TypeError("zone must be a FirstBrillouinZone1D instance.")
        if number_of_points < 2:
            raise ValueError("number_of_points must be at least 2.")
        return cls(
            np.linspace(
                zone.lower,
                zone.upper,
                number_of_points,
                endpoint=endpoint,
                dtype=np.float64,
            )
        )

    @property
    def size(self) -> int:
        """Return the number of sampled wavevectors."""
        return self.values.size

    def check_args(self) -> None:
        """Validate the one-dimensional path."""
        if self.values.ndim != 1 or self.values.size == 0:
            raise ValueError("values must be a nonempty one-dimensional array.")
        if not np.all(np.isfinite(self.values)):
            raise ValueError("values must contain only finite values.")


class KPointPathND:
    """Piecewise-linear path through a multidimensional reciprocal space."""

    def __init__(
        self,
        nodes: FloatArray,
        labels: Sequence[str],
        *,
        points_per_segment: int,
        dimension: int,
    ) -> None:
        self.nodes = np.asarray(nodes, dtype=np.float64)
        self.labels = tuple(labels)
        self.points_per_segment = int(points_per_segment)
        self.dimension = int(dimension)
        self.check_args()
        self.values, self.distances, self.node_distances = self.construct_path()
        self.values.setflags(write=False)
        self.distances.setflags(write=False)
        self.node_distances.setflags(write=False)

    def check_args(self) -> None:
        """Validate path nodes, labels, and sampling density."""
        if self.nodes.ndim != 2 or self.nodes.shape[1] != self.dimension:
            raise ValueError(
                f"nodes must have shape (number_of_nodes, {self.dimension})."
            )
        if self.nodes.shape[0] < 2:
            raise ValueError("at least two path nodes are required.")
        if len(self.labels) != self.nodes.shape[0]:
            raise ValueError("labels must contain one label for each node.")
        if self.points_per_segment < 2:
            raise ValueError("points_per_segment must be at least 2.")
        if not np.all(np.isfinite(self.nodes)):
            raise ValueError("nodes must contain only finite values.")

    def construct_path(self) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Interpolate every segment and construct cumulative path distance."""
        path_parts: list[FloatArray] = []
        distance_parts: list[FloatArray] = []
        node_distances = [0.0]
        cumulative_distance = 0.0

        for segment_index in range(self.nodes.shape[0] - 1):
            start = self.nodes[segment_index]
            stop = self.nodes[segment_index + 1]
            segment = np.linspace(
                start,
                stop,
                self.points_per_segment,
                endpoint=True,
                dtype=np.float64,
            )
            local_distance = np.linalg.norm(segment - start, axis=1)

            # Drop the first point of every segment after the first so shared
            # high-symmetry nodes appear exactly once in the complete path.
            if segment_index > 0:
                segment = segment[1:]
                local_distance = local_distance[1:]

            path_parts.append(segment)
            distance_parts.append(cumulative_distance + local_distance)
            cumulative_distance += float(np.linalg.norm(stop - start))
            node_distances.append(cumulative_distance)

        return (
            np.vstack(path_parts),
            np.concatenate(distance_parts),
            np.asarray(node_distances, dtype=np.float64),
        )

    @property
    def size(self) -> int:
        """Return the number of sampled wavevectors."""
        return self.values.shape[0]


class KPointPath2D(KPointPathND):
    """Piecewise-linear path through a two-dimensional reciprocal space."""

    def __init__(
        self,
        nodes: FloatArray,
        labels: Sequence[str],
        *,
        points_per_segment: int = 51,
    ) -> None:
        super().__init__(
            nodes,
            labels,
            points_per_segment=points_per_segment,
            dimension=2,
        )


class KPointPath3D(KPointPathND):
    """Piecewise-linear path through a three-dimensional reciprocal space."""

    def __init__(
        self,
        nodes: FloatArray,
        labels: Sequence[str],
        *,
        points_per_segment: int = 51,
    ) -> None:
        super().__init__(
            nodes,
            labels,
            points_per_segment=points_per_segment,
            dimension=3,
        )
