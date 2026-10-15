"""Verification of ``HeatMapVisualizer.plot``."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from projectkoios.physkit.visualize.heatmaps import HeatMapVisualizer


def test_maps_zero_to_white_and_sign_to_red_blue() -> None:
    figure, axis = plt.subplots()
    returned = HeatMapVisualizer([[-2.0, 0.0, 2.0]]).plot(
        axis,
        colorbar=False,
    )
    image = axis.images[0]

    assert returned is axis
    np.testing.assert_allclose(
        image.cmap(image.norm(0.0)),
        (1.0, 1.0, 1.0, 1.0),
        rtol=0.0,
        atol=1.0e-12,
    )
    assert image.cmap(image.norm(-2.0))[2] > image.cmap(image.norm(-2.0))[0]
    assert image.cmap(image.norm(2.0))[0] > image.cmap(image.norm(2.0))[2]
    plt.close(figure)


def test_all_zero_matrix_uses_nondegenerate_white_scale() -> None:
    figure, axis = plt.subplots()
    visualizer = HeatMapVisualizer(np.zeros((2, 3)))
    visualizer.plot(axis, colorbar=False)
    image = axis.images[0]

    assert visualizer.color_limit == 1.0
    assert image.norm.vmin == -1.0
    assert image.norm.vmax == 1.0
    np.testing.assert_allclose(
        image.cmap(image.norm(0.0)),
        (1.0, 1.0, 1.0, 1.0),
        rtol=0.0,
        atol=1.0e-12,
    )
    plt.close(figure)


def test_applies_annotations_labels_origin_aspect_and_colorbar() -> None:
    figure, axis = plt.subplots()
    HeatMapVisualizer([[0.0, 1.25], [-0.5, 2.0]]).plot(
        axis,
        annotate=True,
        annotation_format=".2f",
        colorbar=True,
        origin="lower",
        aspect="auto",
        title="matrix",
        xlabel="column",
        ylabel="row",
    )

    assert [text.get_text() for text in axis.texts] == [
        "0.00",
        "1.25",
        "-0.50",
        "2.00",
    ]
    assert axis.images[0].origin == "lower"
    assert axis.get_aspect() == "auto"
    assert axis.get_title() == "matrix"
    assert axis.get_xlabel() == "column"
    assert axis.get_ylabel() == "row"
    assert len(figure.axes) == 2
    plt.close(figure)
