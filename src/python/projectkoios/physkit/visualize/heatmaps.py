"""Generic real-matrix heatmap visualization."""

from __future__ import annotations

from numbers import Real
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from numpy.typing import ArrayLike, NDArray

_WHITE_CENTERED_BLUE_RED = LinearSegmentedColormap.from_list(
    "physkit_white_centered_blue_red",
    ("#0000ff", "#ffffff", "#ff0000"),
    N=257,
)


class HeatMapVisualizer:
    """Visualize a finite real matrix with zero mapped to exact white.

    Color limits are always symmetric about zero. Positive values are red,
    negative values are blue, and zero is white through Matplotlib's ``bwr``
    colormap. An all-zero matrix uses limits ``[-1, 1]`` so its cells remain
    white instead of triggering a degenerate normalization.
    """

    __slots__ = ("_values", "_color_limit")

    def __init__(self, values: ArrayLike) -> None:
        try:
            uncoerced = np.asarray(values, dtype=object)
        except (TypeError, ValueError) as error:
            raise TypeError("values must be a non-ragged real matrix") from error
        if any(isinstance(value, (bool, np.bool_)) for value in uncoerced.flat):
            raise TypeError("values must not contain Boolean entries")

        try:
            array = np.asarray(values)
        except (TypeError, ValueError) as error:
            raise TypeError("values must be a non-ragged real matrix") from error
        if array.dtype.kind not in "iuf":
            raise TypeError("values must contain real numeric entries")
        if array.ndim != 2:
            raise ValueError("values must be two-dimensional")
        if 0 in array.shape:
            raise ValueError("values must have nonzero dimensions")

        canonical = np.array(
            array,
            dtype=np.float64,
            order="C",
            copy=True,
        )
        if not np.all(np.isfinite(canonical)):
            raise ValueError("values must contain only finite entries")
        canonical.setflags(write=False)

        maximum = float(np.max(np.abs(canonical)))
        color_limit = maximum if maximum > 0.0 else 1.0
        self._values: NDArray[np.float64] = canonical
        self._color_limit = color_limit

    @property
    def values(self) -> NDArray[np.float64]:
        """Return an owned C-contiguous copy of the visualized matrix."""
        return np.array(
            self._values,
            dtype=np.float64,
            order="C",
            copy=True,
        )

    @property
    def color_limit(self) -> float:
        """Return the positive symmetric color limit."""
        return self._color_limit

    def plot(
        self,
        ax: Axes | None = None,
        *,
        annotate: bool = False,
        annotation_format: str = "g",
        colorbar: bool = True,
        origin: Literal["upper", "lower"] = "upper",
        aspect: Literal["equal", "auto"] | Real = "equal",
        title: str | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
    ) -> Axes:
        """Plot the matrix and return the containing Matplotlib axes."""
        if not isinstance(annotate, bool):
            raise TypeError("annotate must be a bool")
        if not isinstance(annotation_format, str):
            raise TypeError("annotation_format must be a string")
        if not isinstance(colorbar, bool):
            raise TypeError("colorbar must be a bool")
        if origin not in ("upper", "lower"):
            raise ValueError("origin must be 'upper' or 'lower'")
        if not (
            aspect in ("equal", "auto")
            or (
                not isinstance(aspect, (bool, np.bool_))
                and isinstance(aspect, Real)
                and np.isfinite(float(aspect))
                and float(aspect) > 0.0
            )
        ):
            raise ValueError(
                "aspect must be 'equal', 'auto', or a positive finite real"
            )
        for name, value in (
            ("title", title),
            ("xlabel", xlabel),
            ("ylabel", ylabel),
        ):
            if value is not None and not isinstance(value, str):
                raise TypeError(f"{name} must be a string or None")

        if ax is None:
            _, ax = plt.subplots(figsize=(6.0, 5.0))
        if not isinstance(ax, Axes):
            raise TypeError("ax must be a matplotlib.axes.Axes or None")

        image = ax.imshow(
            self._values,
            cmap=_WHITE_CENTERED_BLUE_RED,
            interpolation="nearest",
            vmin=-self._color_limit,
            vmax=self._color_limit,
            origin=origin,
            aspect=aspect,
        )
        if annotate:
            for (row, column), value in np.ndenumerate(self._values):
                ax.text(
                    column,
                    row,
                    format(float(value), annotation_format),
                    ha="center",
                    va="center",
                    color="black",
                )
        if colorbar:
            ax.figure.colorbar(image, ax=ax, shrink=0.8)
        if title is not None:
            ax.set_title(title)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        return ax
