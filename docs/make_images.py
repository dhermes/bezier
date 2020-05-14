# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Helper to make images that are intended for docs.

To actually execute these functions with the desired inputs, run:

.. code-block:: console

   $ nox -s docs_images
"""

import os

try:
    from matplotlib import patches
    from matplotlib import path as _path_mod
    import matplotlib.pyplot as plt
except ImportError:
    patches = None
    _path_mod = None
    plt = None
import numpy as np

try:
    import seaborn
except ImportError:
    seaborn = None
import bezier
from bezier import _plot_helpers


BLUE = "blue"
GREEN = "green"
RED = "red"
if seaborn is not None:
    seaborn.set()  # Required in ``seaborn >= 0.8``
    # As of ``0.9.0``, this palette has
    # (BLUE, ORANGE, GREEN, RED, PURPLE, BROWN).
    _COLORS = seaborn.color_palette(palette="deep", n_colors=6)
    BLUE = _COLORS[0]
    GREEN = _COLORS[2]
    RED = _COLORS[3]
    del _COLORS
_DOCS_DIR = os.path.abspath(os.path.dirname(__file__))
IMAGES_DIR = os.path.join(_DOCS_DIR, "images")
NO_IMAGES = "GENERATE_IMAGES" not in os.environ


def save_image(figure, filename):
    """Save an image to the docs images directory.

    Args:
        filename (str): The name of the file (not containing
            directory info).
    """
    path = os.path.join(IMAGES_DIR, filename)
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)


def stack1d(*points):
    """Fill out the columns of matrix with a series of points.

    This is because ``np.hstack()`` will just make another 1D vector
    out of them and ``np.vstack()`` will put them in the rows.

    Args:
        points (Tuple[numpy.ndarray, ...]): Tuple of 1D points (i.e.
            arrays with shape ``(2,)``.

    Returns:
        numpy.ndarray: The array with each point in ``points`` as its
        columns.
    """
    result = np.empty((2, len(points)), order="F")
    for index, point in enumerate(points):
        result[:, index] = point
    return result


def linearization_error(nodes):
    """Image for :func:`.linearization_error` docstring."""
    if NO_IMAGES:
        return

    curve = bezier.Curve.from_nodes(nodes)
    line = bezier.Curve.from_nodes(nodes[:, (0, -1)])
    midpoints = np.hstack([curve.evaluate(0.5), line.evaluate(0.5)])
    ax = curve.plot(256, color=BLUE)
    line.plot(256, ax=ax, color=GREEN)
    ax.plot(
        midpoints[0, :], midpoints[1, :], color="black", linestyle="dashed"
    )
    ax.axis("scaled")
    save_image(ax.figure, "linearization_error.png")


def newton_refine1(s, new_s, curve1, t, new_t, curve2):
    """Image for :func:`.newton_refine` docstring."""
    if NO_IMAGES:
        return

    points = np.hstack([curve1.evaluate(s), curve2.evaluate(t)])
    points_new = np.hstack([curve1.evaluate(new_s), curve2.evaluate(new_t)])
    ax = curve1.plot(256, color=BLUE)
    curve2.plot(256, ax=ax, color=GREEN)
    ax.plot(
        points[0, :],
        points[1, :],
        color="black",
        linestyle="None",
        marker="o",
        markeredgewidth=1,
        markerfacecolor="None",
    )
    ax.plot(
        points_new[0, :],
        points_new[1, :],
        color="black",
        linestyle="None",
        marker="o",
    )
    ax.axis("scaled")
    save_image(ax.figure, "newton_refine1.png")


def newton_refine2(s_vals, curve1, curve2):
    """Image for :func:`.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve1.plot(256, color=BLUE)
    ax.lines[-1].zorder = 1
    curve2.plot(256, ax=ax, color=GREEN)
    ax.lines[-1].zorder = 1
    points = curve1.evaluate_multi(np.asfortranarray(s_vals))
    colors = seaborn.dark_palette("blue", 5)
    ax.scatter(
        points[0, :], points[1, :], c=colors, s=20, alpha=0.75, zorder=2
    )
    ax.axis("scaled")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    save_image(ax.figure, "newton_refine2.png")


def newton_refine3(s_vals, curve1, curve2):
    """Image for :func:`.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve1.plot(256, color=BLUE)
    ax.lines[-1].zorder = 1
    curve2.plot(256, ax=ax, color=GREEN)
    ax.lines[-1].zorder = 1
    points = curve1.evaluate_multi(np.asfortranarray(s_vals))
    colors = seaborn.dark_palette("blue", 6)
    ax.scatter(
        points[0, :], points[1, :], c=colors, s=20, alpha=0.75, zorder=2
    )
    ax.axis("scaled")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.5625)
    save_image(ax.figure, "newton_refine3.png")


def segment_intersection1(start0, end0, start1, end1, s):
    """Image for :func:`.segment_intersection` docstring."""
    if NO_IMAGES:
        return

    line0 = bezier.Curve.from_nodes(stack1d(start0, end0))
    line1 = bezier.Curve.from_nodes(stack1d(start1, end1))
    ax = line0.plot(2, color=BLUE)
    line1.plot(256, ax=ax, color=GREEN)
    (x_val,), (y_val,) = line0.evaluate(s)
    ax.plot([x_val], [y_val], color="black", marker="o")
    ax.axis("scaled")
    save_image(ax.figure, "segment_intersection1.png")


def segment_intersection2(start0, end0, start1, end1):
    """Image for :func:`.segment_intersection` docstring."""
    if NO_IMAGES:
        return

    line0 = bezier.Curve.from_nodes(stack1d(start0, end0))
    line1 = bezier.Curve.from_nodes(stack1d(start1, end1))
    ax = line0.plot(2, color=BLUE)
    line1.plot(2, ax=ax, color=GREEN)
    ax.axis("scaled")
    save_image(ax.figure, "segment_intersection2.png")


def helper_parallel_lines(start0, end0, start1, end1, filename):
    """Image for :func:`.parallel_lines_parameters` docstring."""
    if NO_IMAGES:
        return

    figure = plt.figure()
    ax = figure.gca()
    points = stack1d(start0, end0, start1, end1)
    ax.plot(points[0, :2], points[1, :2], marker="o", color=BLUE)
    ax.plot(points[0, 2:], points[1, 2:], marker="o", color=GREEN)
    ax.axis("scaled")
    _plot_helpers.add_plot_boundary(ax)
    save_image(figure, filename)


def add_patch(
    ax, nodes, color, with_nodes=True, alpha=0.625, node_color="black"
):
    # ``nodes`` is stored Fortran-contiguous with ``x-y`` points in each
    # column but ``Path()`` wants ``x-y`` points in each row.
    path = _path_mod.Path(nodes.T)
    patch = patches.PathPatch(
        path, edgecolor=color, facecolor=color, alpha=alpha
    )
    ax.add_patch(patch)
    if with_nodes:
        ax.plot(
            nodes[0, :],
            nodes[1, :],
            color=node_color,
            linestyle="None",
            marker="o",
        )


def curve_constructor(curve):
    """Image for :class`.Curve` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    line = ax.lines[0]
    nodes = curve._nodes
    ax.plot(
        nodes[0, :], nodes[1, :], color="black", linestyle="None", marker="o"
    )
    add_patch(ax, nodes, line.get_color())
    ax.axis("scaled")
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, "curve_constructor.png")


def curve_evaluate(curve):
    """Image for :meth`.Curve.evaluate` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    points = curve.evaluate_multi(np.asfortranarray([0.75]))
    ax.plot(
        points[0, :], points[1, :], color="black", linestyle="None", marker="o"
    )
    ax.axis("scaled")
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, "curve_evaluate.png")


def curve_evaluate_hodograph(curve, s):
    """Image for :meth`.Curve.evaluate_hodograph` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)

    points = curve.evaluate_multi(np.asfortranarray([s]))
    if points.shape != (2, 1):
        raise ValueError("Unexpected shape", points)
    point = points[:, 0]

    tangents = curve.evaluate_hodograph(s)
    if tangents.shape != (2, 1):
        raise ValueError("Unexpected shape", tangents)
    tangent = tangents[:, 0]

    ax.plot(
        [point[0] - 2 * tangent[0], point[0] + 2 * tangent[0]],
        [point[1] - 2 * tangent[1], point[1] + 2 * tangent[1]],
        color=BLUE,
        alpha=0.5,
    )

    ax.plot(
        [point[0], point[0] + tangent[0]],
        [point[1], point[1] + tangent[1]],
        color="black",
        linestyle="dashed",
        marker="o",
        markersize=5,
    )

    ax.axis("scaled")
    ax.set_xlim(-0.125, 1.75)
    ax.set_ylim(-0.0625, 0.75)
    save_image(ax.figure, "curve_evaluate_hodograph.png")


def curve_subdivide(curve, left, right):
    """Image for :meth`.Curve.subdivide` docstring."""
    if NO_IMAGES:
        return

    figure = plt.figure()
    ax = figure.gca()
    add_patch(ax, curve._nodes, "gray")
    ax = left.plot(256, ax=ax, color=BLUE)
    line = ax.lines[-1]
    add_patch(ax, left._nodes, line.get_color())
    right.plot(256, ax=ax, color=GREEN)
    line = ax.lines[-1]
    add_patch(ax, right._nodes, line.get_color())
    ax.axis("scaled")
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 3.125)
    save_image(ax.figure, "curve_subdivide.png")


def curve_intersect(curve1, curve2, s_vals):
    """Image for :meth`.Curve.intersect` docstring."""
    if NO_IMAGES:
        return

    ax = curve1.plot(256, color=BLUE)
    curve2.plot(256, ax=ax, color=GREEN)
    intersections = curve1.evaluate_multi(s_vals)
    ax.plot(
        intersections[0, :],
        intersections[1, :],
        color="black",
        linestyle="None",
        marker="o",
    )
    ax.axis("scaled")
    ax.set_xlim(0.0, 0.75)
    ax.set_ylim(0.0, 0.75)
    save_image(ax.figure, "curve_intersect.png")


def triangle_constructor(triangle):
    """Image for :class`.Triangle` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE, with_nodes=True)
    line = ax.lines[0]
    nodes = triangle._nodes
    add_patch(ax, nodes[:, (0, 1, 2, 5)], line.get_color())
    delta = 1.0 / 32.0
    ax.text(
        nodes[0, 0],
        nodes[1, 0],
        r"$v_0$",
        fontsize=20,
        verticalalignment="top",
        horizontalalignment="right",
    )
    ax.text(
        nodes[0, 1],
        nodes[1, 1],
        r"$v_1$",
        fontsize=20,
        verticalalignment="top",
        horizontalalignment="center",
    )
    ax.text(
        nodes[0, 2],
        nodes[1, 2],
        r"$v_2$",
        fontsize=20,
        verticalalignment="top",
        horizontalalignment="left",
    )
    ax.text(
        nodes[0, 3] - delta,
        nodes[1, 3],
        r"$v_3$",
        fontsize=20,
        verticalalignment="center",
        horizontalalignment="right",
    )
    ax.text(
        nodes[0, 4] + delta,
        nodes[1, 4],
        r"$v_4$",
        fontsize=20,
        verticalalignment="center",
        horizontalalignment="left",
    )
    ax.text(
        nodes[0, 5],
        nodes[1, 5] + delta,
        r"$v_5$",
        fontsize=20,
        verticalalignment="bottom",
        horizontalalignment="center",
    )
    ax.axis("scaled")
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, "triangle_constructor.png")


def triangle_evaluate_barycentric(triangle, point):
    """Image for :meth`.Triangle.evaluate_barycentric` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE)
    ax.plot(
        point[0, :], point[1, :], color="black", linestyle="None", marker="o"
    )
    ax.axis("scaled")
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, "triangle_evaluate_barycentric.png")


def triangle_evaluate_cartesian_multi(triangle, points):
    """Image for :meth`.Triangle.evaluate_cartesian_multi` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE)
    ax.plot(
        points[0, :], points[1, :], color="black", linestyle="None", marker="o"
    )
    delta = 1.0 / 32.0
    font_size = 18
    ax.text(
        points[0, 0],
        points[1, 0],
        r"$w_0$",
        fontsize=font_size,
        verticalalignment="top",
        horizontalalignment="right",
    )
    ax.text(
        points[0, 1] + 2 * delta,
        points[1, 1],
        r"$w_1$",
        fontsize=font_size,
        verticalalignment="center",
        horizontalalignment="left",
    )
    ax.text(
        points[0, 2],
        points[1, 2] + delta,
        r"$w_2$",
        fontsize=font_size,
        verticalalignment="bottom",
        horizontalalignment="left",
    )
    ax.axis("scaled")
    ax.set_xlim(-3.125, 2.375)
    ax.set_ylim(-0.25, 2.125)
    save_image(ax.figure, "triangle_evaluate_cartesian_multi.png")


def triangle_evaluate_barycentric_multi(triangle, points):
    """Image for :meth`.Triangle.evaluate_barycentric_multi` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE)
    ax.plot(
        points[0, :], points[1, :], color="black", linestyle="None", marker="o"
    )
    delta = 1.0 / 32.0
    font_size = 18
    ax.text(
        points[0, 0],
        points[1, 0] + delta,
        r"$w_0$",
        fontsize=font_size,
        verticalalignment="bottom",
        horizontalalignment="center",
    )
    ax.text(
        points[0, 1],
        points[1, 1] - delta,
        r"$w_1$",
        fontsize=font_size,
        verticalalignment="top",
        horizontalalignment="right",
    )
    ax.text(
        points[0, 2],
        points[1, 2],
        r"$w_2$",
        fontsize=font_size,
        verticalalignment="bottom",
        horizontalalignment="left",
    )
    ax.text(
        points[0, 3],
        points[1, 3],
        r"$w_3$",
        fontsize=font_size,
        verticalalignment="top",
        horizontalalignment="right",
    )
    ax.axis("scaled")
    ax.set_xlim(-3.125, 2.125)
    ax.set_ylim(-0.3125, 2.125)
    save_image(ax.figure, "triangle_evaluate_barycentric_multi.png")


def triangle_is_valid1(triangle):
    """Image for :meth`.Triangle.is_valid` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE)
    ax.axis("scaled")
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 2.125)
    save_image(ax.figure, "triangle_is_valid1.png")


def triangle_is_valid2(triangle):
    """Image for :meth`.Triangle.is_valid` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE)
    ax.axis("scaled")
    ax.set_xlim(-0.125, 1.0625)
    ax.set_ylim(-0.0625, 1.0625)
    save_image(ax.figure, "triangle_is_valid2.png")


def triangle_is_valid3(triangle):
    """Image for :meth`.Triangle.is_valid` docstring."""
    if NO_IMAGES:
        return

    edge1, edge2, edge3 = triangle.edges
    N = 128
    # Compute points on each edge.
    std_s = np.linspace(0.0, 1.0, N + 1)
    points1 = edge1.evaluate_multi(std_s)
    points2 = edge2.evaluate_multi(std_s)
    points3 = edge3.evaluate_multi(std_s)
    # Compute the actual boundary where the Jacobian is 0.
    s_vals = np.linspace(0.0, 0.2, N)
    t_discrim = np.sqrt((1.0 - s_vals) * (1.0 - 5.0 * s_vals))
    t_top = 0.5 * (1.0 - s_vals + t_discrim)
    t_bottom = 0.5 * (1.0 - s_vals - t_discrim)
    jacobian_zero_params = np.zeros((2 * N - 1, 2), order="F")
    jacobian_zero_params[:N, 0] = s_vals
    jacobian_zero_params[:N, 1] = t_top
    jacobian_zero_params[N:, 0] = s_vals[-2::-1]
    jacobian_zero_params[N:, 1] = t_bottom[-2::-1]
    jac_edge = triangle.evaluate_cartesian_multi(jacobian_zero_params)
    # Add the triangle to the plot and add a dashed line
    # for each "true" edge.
    figure = plt.figure()
    ax = figure.gca()
    (line,) = ax.plot(jac_edge[0, :], jac_edge[1, :], color=BLUE)
    color = line.get_color()
    ax.plot(points1[0, :], points1[1, :], color="black", linestyle="dashed")
    ax.plot(points2[0, :], points2[1, :], color="black", linestyle="dashed")
    ax.plot(points3[0, :], points3[1, :], color="black", linestyle="dashed")
    polygon = np.hstack([points1[:, 1:], points2[:, 1:], jac_edge[:, 1:]])
    add_patch(ax, polygon, color, with_nodes=False)
    ax.axis("scaled")
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.0625, 1.0625)
    save_image(ax.figure, "triangle_is_valid3.png")


def triangle_subdivide1():
    """Image for :meth`.Triangle.subdivide` docstring."""
    if NO_IMAGES:
        return

    triangle = bezier.Triangle.from_nodes(
        np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    )
    triangle_a, triangle_b, triangle_c, triangle_d = triangle.subdivide()
    figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    for ax in (ax1, ax2, ax3, ax4):
        triangle.plot(2, ax=ax, color=BLUE)
    triangle_a.plot(2, ax=ax1, color=GREEN)
    ax1.text(
        1.0 / 6.0,
        1.0 / 6.0,
        r"$A$",
        fontsize=20,
        verticalalignment="center",
        horizontalalignment="center",
    )
    triangle_b.plot(2, ax=ax2, color=GREEN)
    ax2.text(
        1.0 / 3.0,
        1.0 / 3.0,
        r"$B$",
        fontsize=20,
        verticalalignment="center",
        horizontalalignment="center",
    )
    triangle_c.plot(2, ax=ax3, color=GREEN)
    ax3.text(
        2.0 / 3.0,
        1.0 / 6.0,
        r"$C$",
        fontsize=20,
        verticalalignment="center",
        horizontalalignment="center",
    )
    triangle_d.plot(2, ax=ax4, color=GREEN)
    ax4.text(
        1.0 / 6.0,
        2.0 / 3.0,
        r"$D$",
        fontsize=20,
        verticalalignment="center",
        horizontalalignment="center",
    )
    for ax in (ax1, ax2, ax3, ax4):
        ax.axis("scaled")
    save_image(figure, "triangle_subdivide1")


def add_edges(ax, triangle, s_vals, color):
    edge1, edge2, edge3 = triangle.edges
    # Compute points on each edge.
    points1 = edge1.evaluate_multi(s_vals)
    points2 = edge2.evaluate_multi(s_vals)
    points3 = edge3.evaluate_multi(s_vals)
    # Add the points to the plot.
    ax.plot(points1[0, :], points1[1, :], color=color)
    ax.plot(points2[0, :], points2[1, :], color=color)
    ax.plot(points3[0, :], points3[1, :], color=color)


def triangle_subdivide2(triangle, sub_triangle_b):
    """Image for :meth`.Triangle.subdivide` docstring."""
    if NO_IMAGES:
        return

    # Plot set-up.
    figure = plt.figure()
    ax = figure.gca()
    colors = seaborn.husl_palette(6)
    N = 128
    s_vals = np.linspace(0.0, 1.0, N + 1)
    # Add edges from triangle.
    add_edges(ax, triangle, s_vals, colors[4])
    # Now do the same for triangle B.
    add_edges(ax, sub_triangle_b, s_vals, colors[0])
    # Add the control points polygon for the original triangle.
    nodes = triangle._nodes[:, (0, 2, 4, 5, 0)]
    add_patch(ax, nodes, colors[2], with_nodes=False)
    # Add the control points polygon for the sub-triangle.
    nodes = sub_triangle_b._nodes[:, (0, 1, 2, 5, 3, 0)]
    add_patch(ax, nodes, colors[1], with_nodes=False)
    # Plot **all** the nodes.
    sub_nodes = sub_triangle_b._nodes
    ax.plot(
        sub_nodes[0, :],
        sub_nodes[1, :],
        color="black",
        linestyle="None",
        marker="o",
    )
    # Take those same points and add the boundary.
    ax.plot(nodes[0, :], nodes[1, :], color="black", linestyle="dashed")
    ax.axis("scaled")
    ax.set_xlim(-1.125, 2.125)
    ax.set_ylim(-0.125, 4.125)
    save_image(ax.figure, "triangle_subdivide2")


def curved_polygon_constructor1(curved_poly):
    """Image for :class`.CurvedPolygon` docstring."""
    if NO_IMAGES:
        return

    ax = curved_poly.plot(256, color=BLUE)
    ax.axis("scaled")
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.625, 1.625)
    save_image(ax.figure, "curved_polygon_constructor1.png")


def curved_polygon_constructor2(curved_poly):
    """Image for :class`.CurvedPolygon` docstring."""
    if NO_IMAGES:
        return

    ax = curved_poly.plot(256, color=BLUE)
    ax.axis("scaled")
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, "curved_polygon_constructor2.png")


def triangle_locate(triangle, point):
    """Image for :meth`.Triangle.locate` docstring."""
    if NO_IMAGES:
        return

    ax = triangle.plot(256, color=BLUE)
    ax.plot(
        point[0, :], point[1, :], color="black", linestyle="None", marker="o"
    )
    ax.axis("scaled")
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.1875, 1.0625)
    save_image(ax.figure, "triangle_locate.png")


def curve_specialize(curve, new_curve):
    """Image for :meth`.Curve.specialize` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    interval = r"$\left[0, 1\right]$"
    line = ax.lines[-1]
    line.set_label(interval)
    color1 = line.get_color()
    new_curve.plot(256, ax=ax, color=GREEN)
    interval = r"$\left[-\frac{1}{4}, \frac{3}{4}\right]$"
    line = ax.lines[-1]
    line.set_label(interval)
    ax.plot(
        curve._nodes[0, (0, -1)],
        curve._nodes[1, (0, -1)],
        color=color1,
        linestyle="None",
        marker="o",
    )
    ax.plot(
        new_curve._nodes[0, (0, -1)],
        new_curve._nodes[1, (0, -1)],
        color=line.get_color(),
        linestyle="None",
        marker="o",
    )
    ax.legend(loc="lower right", fontsize=12)
    ax.axis("scaled")
    ax.set_xlim(-0.375, 1.125)
    ax.set_ylim(-0.75, 0.625)
    save_image(ax.figure, "curve_specialize.png")


def newton_refine_triangle(triangle, x_val, y_val, s, t, new_s, new_t):
    """Image for :func:`.hazmat.triangle_helpers.newton_refine` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2)
    # Plot features of the parameter space in ax1.
    linear_triangle = bezier.Triangle.from_nodes(
        np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    )
    linear_triangle.plot(2, ax=ax1, color=BLUE)
    ax1.plot([0.25], [0.5], marker="H", color=GREEN)
    ax1.plot([s], [t], color="black", linestyle="None", marker="o")
    ax1.plot(
        [new_s],
        [new_t],
        color="black",
        linestyle="None",
        marker="o",
        markeredgewidth=1,
        markerfacecolor="None",
    )
    # Plot the equivalent output in ax2.
    triangle.plot(256, ax=ax2, color=BLUE)
    points = triangle.evaluate_cartesian_multi(
        np.asfortranarray([[s, t], [new_s, new_t]])
    )
    ax2.plot([x_val], [y_val], marker="H", color=GREEN)
    ax2.plot(
        points[0, [0]],
        points[1, [0]],
        color="black",
        linestyle="None",
        marker="o",
    )
    ax2.plot(
        points[0, [1]],
        points[1, [1]],
        color="black",
        linestyle="None",
        marker="o",
        markeredgewidth=1,
        markerfacecolor="None",
    )
    # Set the axis bounds / scaling.
    ax1.axis("scaled")
    ax1.set_xlim(-0.0625, 1.0625)
    ax1.set_ylim(-0.0625, 1.0625)
    ax2.axis("scaled")
    ax2.set_xlim(-0.125, 2.125)
    ax2.set_ylim(-0.125, 2.125)
    save_image(figure, "newton_refine_triangle.png")


def classify_help(s, curve1, triangle1, curve2, triangle2, interior, ax=None):
    assert triangle1.is_valid
    edge1, _, _ = triangle1.edges
    assert np.all(edge1._nodes == curve1._nodes)
    assert triangle2.is_valid
    edge2, _, _ = triangle2.edges
    assert np.all(edge2._nodes == curve2._nodes)
    ax = triangle1.plot(256, ax=ax, color=BLUE)
    # Manually reduce the alpha on the triangle patch(es).
    ax.patches[-1].set_alpha(0.1875)
    color1 = ax.lines[-1].get_color()
    triangle2.plot(256, ax=ax, color=GREEN)
    ax.patches[-1].set_alpha(0.1875)
    color2 = ax.lines[-1].get_color()
    # Remove the existing boundary (lines) and just add our edges.
    while ax.lines:
        ax.lines[-1].remove()
    edge1.plot(256, ax=ax, color=color1)
    edge2.plot(256, ax=ax, color=color2)
    (int_x,), (int_y,) = curve1.evaluate(s)
    if interior == 0:
        color = color1
    elif interior == 1:
        color = color2
    else:
        color = RED
    ax.plot([int_x], [int_y], color=color, linestyle="None", marker="o")
    ax.axis("scaled")
    return ax


def classify_intersection1(s, curve1, tangent1, curve2, tangent2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [[1.0, 1.75, 2.0, 1.0, 1.5, 1.0], [0.0, 0.25, 1.0, 1.0, 1.5, 2.0]]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [0.0, 1.6875, 2.0, 0.25, 1.25, 0.5],
                [0.0, 0.0625, 0.5, 1.0, 1.25, 2.0],
            ]
        )
    )
    ax = classify_help(s, curve1, triangle1, curve2, triangle2, 0)
    (int_x,), (int_y,) = curve1.evaluate(s)
    # Remove the alpha from the color
    color1 = ax.patches[0].get_facecolor()[:3]
    color2 = ax.patches[1].get_facecolor()[:3]
    ax.plot(
        [int_x, int_x + tangent1[0, 0]],
        [int_y, int_y + tangent1[1, 0]],
        color=color1,
        linestyle="dashed",
    )
    ax.plot(
        [int_x, int_x + tangent2[0, 0]],
        [int_y, int_y + tangent2[1, 0]],
        color=color2,
        linestyle="dashed",
    )
    ax.plot([int_x], [int_y], color=color1, linestyle="None", marker="o")
    ax.axis("scaled")
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, "classify_intersection1.png")


def classify_intersection2(s, curve1, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [[1.0, 1.5, 2.0, 1.25, 1.75, 1.5], [0.0, 1.0, 0.0, 1.0, 1.0, 2.0]]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [[0.0, 1.5, 3.0, 0.75, 2.25, 1.5], [0.0, 1.0, 0.0, 2.0, 2.0, 4.0]]
        )
    )
    ax = classify_help(s, curve1, triangle1, curve2, triangle2, 1)
    ax.set_xlim(-0.0625, 3.0625)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, "classify_intersection2.png")


def classify_intersection3(s, curve1, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [2.0, 1.5, 1.0, 1.75, 1.25, 1.5],
                [0.0, 1.0, 0.0, -1.0, -1.0, -2.0],
            ]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [3.0, 1.5, 0.0, 2.25, 0.75, 1.5],
                [0.0, 1.0, 0.0, -2.0, -2.0, -4.0],
            ]
        )
    )
    ax = classify_help(s, curve1, triangle1, curve2, triangle2, 0)
    ax.set_xlim(-0.0625, 3.0625)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, "classify_intersection3.png")


def classify_intersection4(s, curve1, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [2.0, 1.5, 1.0, 1.75, 1.25, 1.5],
                [0.0, 1.0, 0.0, -1.0, -1.0, -2.0],
            ]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [[0.0, 1.5, 3.0, 0.75, 2.25, 1.5], [0.0, 1.0, 0.0, 2.0, 2.0, 4.0]]
        )
    )
    ax = classify_help(s, curve1, triangle1, curve2, triangle2, None)
    ax.set_xlim(-0.0625, 3.0625)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, "classify_intersection4.png")


def classify_intersection5(s, curve1, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [1.0, 1.5, 2.0, 1.25, 1.75, 1.5],
                [0.0, 1.0, 0.0, 0.9375, 0.9375, 1.875],
            ]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [3.0, 1.5, 0.0, 2.25, 0.75, 1.5],
                [0.0, 1.0, 0.0, -2.0, -2.0, -4.0],
            ]
        )
    )
    figure, (ax1, ax2) = plt.subplots(2, 1)
    classify_help(s, curve1, triangle1, curve2, triangle2, 0, ax=ax1)
    classify_help(s, curve1, triangle1, curve2, triangle2, 1, ax=ax2)
    # Remove the alpha from the color
    color1 = ax1.patches[0].get_facecolor()[:3]
    color2 = ax1.patches[1].get_facecolor()[:3]
    # Now add the "degenerate" intersection polygons. The first
    # comes from specializing to
    # left1(0.5, 1.0)-left2(0.0, 0.25)-right1(0.375, 0.5)
    triangle3 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [1.5, 1.75, 2.0, 1.6875, 1.9375, 1.875],
                [0.5, 0.5, 0.0, 0.5, 0.234375, 0.46875],
            ]
        )
    )
    # NOTE: We don't require the intersection polygon be valid.
    triangle3.plot(256, ax=ax1, color=RED)
    # The second comes from specializing to
    # left1(0.0, 0.5)-right1(0.5, 0.625)-left3(0.75, 1.0)
    triangle4 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [1.0, 1.25, 1.5, 1.0625, 1.3125, 1.125],
                [0.0, 0.5, 0.5, 0.234375, 0.5, 0.46875],
            ]
        )
    )
    # NOTE: We don't require the intersection polygon be valid.
    triangle4.plot(256, ax=ax2, color=RED)
    (int_x,), (int_y,) = curve1.evaluate(s)
    ax1.plot([int_x], [int_y], color=color1, linestyle="None", marker="o")
    ax2.plot([int_x], [int_y], color=color2, linestyle="None", marker="o")
    for ax in (ax1, ax2):
        ax.axis("scaled")
        ax.set_xlim(-0.0625, 3.0625)
        ax.set_ylim(-0.0625, 0.5625)
    plt.setp(ax1.get_xticklabels(), visible=False)
    figure.tight_layout(h_pad=-7.0)
    save_image(figure, "classify_intersection5.png")


def classify_intersection6(s, curve1, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [-0.125, -0.125, 0.375, -0.0625, 0.1875, 0.0],
                [0.0625, -0.0625, 0.0625, 0.15625, 0.15625, 0.25],
            ]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [-0.25, -0.25, 0.75, 0.125, 0.625, 0.5],
                [0.25, -0.25, 0.25, 0.625, 0.625, 1.0],
            ]
        )
    )
    ax = classify_help(s, curve1, triangle1, curve2, triangle2, None)
    ax.set_xlim(-0.3125, 1.0625)
    ax.set_ylim(-0.0625, 0.3125)
    save_image(ax.figure, "classify_intersection6.png")


def classify_intersection7(s, curve1a, curve1b, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [0.0, 4.5, 9.0, 0.0, 4.5, 0.0],
                [0.0, 0.0, 2.25, 1.25, 2.375, 2.5],
            ]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [11.25, 9.0, 2.75, 8.125, 3.875, 5.0],
                [0.0, 4.5, 1.0, -0.75, -0.25, -1.5],
            ]
        )
    )
    figure, (ax1, ax2) = plt.subplots(2, 1)
    classify_help(s, curve1a, triangle1, curve2, triangle2, None, ax=ax1)
    triangle1._nodes = np.asfortranarray(
        triangle1._nodes[:, (2, 4, 5, 1, 3, 0)]
    )
    triangle1._edges = None
    classify_help(0.0, curve1b, triangle1, curve2, triangle2, 0, ax=ax2)
    for ax in (ax1, ax2):
        ax.set_xlim(-0.125, 11.5)
        ax.set_ylim(-0.125, 2.625)
    plt.setp(ax1.get_xticklabels(), visible=False)
    figure.tight_layout(h_pad=-5.0)
    save_image(figure, "classify_intersection7.png")


def get_curvature(nodes, s, tangent_vec, curvature):
    """Image for :func:`get_curvature` docstring."""
    if NO_IMAGES:
        return

    curve = bezier.Curve.from_nodes(nodes)
    # Find the center of the circle along the direction
    # perpendicular to the tangent vector (90 degree left turn).
    radius_dir = np.asfortranarray([[-tangent_vec[1, 0]], [tangent_vec[0, 0]]])
    radius_dir /= np.linalg.norm(radius_dir, ord=2)
    point = curve.evaluate(s)
    circle_center = point + radius_dir / curvature
    # Add the curve.
    ax = curve.plot(256, color=BLUE)
    # Add the circle.
    circle_center = circle_center.ravel(order="F")
    circle = plt.Circle(circle_center, 1.0 / abs(curvature), alpha=0.25)
    ax.add_artist(circle)
    # Add the point.
    ax.plot(
        point[0, :], point[1, :], color="black", marker="o", linestyle="None"
    )
    ax.axis("scaled")
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.0625, 0.625)
    save_image(ax.figure, "get_curvature.png")


def curve_locate(curve, point1, point2, point3):
    """Image for :meth`.Curve.locate` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    points = np.hstack([point1, point2, point3])
    ax.plot(
        points[0, :], points[1, :], color="black", linestyle="None", marker="o"
    )
    ax.axis("scaled")
    ax.set_xlim(-0.8125, 0.0625)
    ax.set_ylim(0.75, 2.0625)
    save_image(ax.figure, "curve_locate.png")


def newton_refine_curve(curve, point, s, new_s):
    """Image for :func:`.hazmat.curve_helpers.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    ax.plot(point[0, :], point[1, :], marker="H", color=GREEN)
    wrong_points = curve.evaluate_multi(np.asfortranarray([s, new_s]))
    ax.plot(
        wrong_points[0, [0]],
        wrong_points[1, [0]],
        color="black",
        linestyle="None",
        marker="o",
    )
    ax.plot(
        wrong_points[0, [1]],
        wrong_points[1, [1]],
        color="black",
        linestyle="None",
        marker="o",
        markeredgewidth=1,
        markerfacecolor="None",
    )
    # Set the axis bounds / scaling.
    ax.axis("scaled")
    ax.set_xlim(-0.125, 3.125)
    ax.set_ylim(-0.125, 1.375)
    save_image(ax.figure, "newton_refine_curve.png")


def newton_refine_curve_cusp(curve, s_vals):
    """Image for :func:`.hazmat.curve_helpers.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    ax.lines[-1].zorder = 1
    points = curve.evaluate_multi(np.asfortranarray(s_vals))
    colors = seaborn.dark_palette("blue", 6)
    ax.scatter(
        points[0, :], points[1, :], c=colors, s=20, alpha=0.75, zorder=2
    )
    # Set the axis bounds / scaling.
    ax.axis("scaled")
    ax.set_xlim(-0.125, 6.125)
    ax.set_ylim(-3.125, 3.125)
    save_image(ax.figure, "newton_refine_curve_cusp.png")


def classify_intersection8(s, curve1, triangle1, curve2, triangle2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    ax = classify_help(s, curve1, triangle1, curve2, triangle2, None)
    ax.set_xlim(-1.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, "classify_intersection8.png")


def _edges_classify_intersection9():
    """The edges for the curved polygon intersection used below.

    Helper for :func:`classify_intersection9`.
    """
    edges1 = (
        bezier.Curve.from_nodes(
            np.asfortranarray([[32.0, 30.0], [20.0, 25.0]])
        ),
        bezier.Curve.from_nodes(
            np.asfortranarray([[30.0, 25.0, 20.0], [25.0, 20.0, 20.0]])
        ),
        bezier.Curve.from_nodes(
            np.asfortranarray([[20.0, 25.0, 30.0], [20.0, 20.0, 15.0]])
        ),
        bezier.Curve.from_nodes(
            np.asfortranarray([[30.0, 32.0], [15.0, 20.0]])
        ),
    )
    edges2 = (
        bezier.Curve.from_nodes(
            np.asfortranarray([[8.0, 10.0], [20.0, 15.0]])
        ),
        bezier.Curve.from_nodes(
            np.asfortranarray([[10.0, 15.0, 20.0], [15.0, 20.0, 20.0]])
        ),
        bezier.Curve.from_nodes(
            np.asfortranarray([[20.0, 15.0, 10.0], [20.0, 20.0, 25.0]])
        ),
        bezier.Curve.from_nodes(
            np.asfortranarray([[10.0, 8.0], [25.0, 20.0]])
        ),
    )
    return edges1, edges2


def classify_intersection9(s, curve1, curve2):
    """Image for :func:`.hazmat.triangle_helpers.classify_intersection` doc."""
    if NO_IMAGES:
        return

    triangle1 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [0.0, 20.0, 40.0, 10.0, 30.0, 20.0],
                [0.0, 40.0, 0.0, 25.0, 25.0, 50.0],
            ]
        )
    )
    triangle2 = bezier.Triangle.from_nodes(
        np.asfortranarray(
            [
                [40.0, 20.0, 0.0, 30.0, 10.0, 20.0],
                [40.0, 0.0, 40.0, 15.0, 15.0, -10.0],
            ]
        )
    )
    figure, (ax1, ax2) = plt.subplots(1, 2)
    classify_help(s, curve1, triangle1, curve2, triangle2, 0, ax=ax1)
    classify_help(s, curve1, triangle1, curve2, triangle2, 1, ax=ax2)
    # Remove the alpha from the color
    color1 = ax1.patches[0].get_facecolor()[:3]
    color2 = ax1.patches[1].get_facecolor()[:3]
    # Now add the "degenerate" intersection polygons.
    cp_edges1, cp_edges2 = _edges_classify_intersection9()
    curved_polygon1 = bezier.CurvedPolygon(*cp_edges1)
    curved_polygon1.plot(256, ax=ax1, color=RED)
    curved_polygon2 = bezier.CurvedPolygon(*cp_edges2)
    curved_polygon2.plot(256, ax=ax2, color=RED)
    (int_x,), (int_y,) = curve1.evaluate(s)
    ax1.plot([int_x], [int_y], color=color1, linestyle="None", marker="o")
    ax2.plot([int_x], [int_y], color=color2, linestyle="None", marker="o")
    for ax in (ax1, ax2):
        ax.axis("scaled")
        ax.set_xlim(-2.0, 42.0)
        ax.set_ylim(-12.0, 52.0)
    plt.setp(ax2.get_yticklabels(), visible=False)
    figure.tight_layout(w_pad=1.0)
    save_image(figure, "classify_intersection9.png")


def curve_elevate(curve, elevated):
    """Image for :meth:`.curve.Curve.elevate` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2)
    curve.plot(256, ax=ax1, color=BLUE)
    color = ax1.lines[-1].get_color()
    add_patch(ax1, curve._nodes, color)
    elevated.plot(256, ax=ax2, color=BLUE)
    color = ax2.lines[-1].get_color()
    add_patch(ax2, elevated._nodes, color)
    ax1.axis("scaled")
    ax2.axis("scaled")
    _plot_helpers.add_plot_boundary(ax1)
    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*ax1.get_ylim())
    save_image(figure, "curve_elevate.png")


def triangle_elevate(triangle, elevated):
    """Image for :meth:`.triangle.Triangle.elevate` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2)
    triangle.plot(256, ax=ax1, color=BLUE)
    color = ax1.lines[-1].get_color()
    nodes = triangle._nodes[:, (0, 1, 2, 4, 5)]
    add_patch(ax1, nodes, color)
    elevated.plot(256, ax=ax2, color=BLUE)
    color = ax2.lines[-1].get_color()
    nodes = elevated._nodes[:, (0, 1, 2, 3, 6, 8, 9)]
    add_patch(ax2, nodes, color)
    ax1.axis("scaled")
    ax2.axis("scaled")
    _plot_helpers.add_plot_boundary(ax1)
    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*ax1.get_ylim())
    save_image(figure, "triangle_elevate.png")


def unit_triangle():
    """Image for :class:`.triangle.Triangle` docstring."""
    if NO_IMAGES:
        return

    nodes = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    triangle = bezier.Triangle(nodes, degree=1)
    ax = triangle.plot(256, color=BLUE)
    ax.axis("scaled")
    _plot_helpers.add_plot_boundary(ax)
    save_image(ax.figure, "unit_triangle.png")


def curve_reduce(curve, reduced):
    """Image for :meth:`.curve.Curve.reduce` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
    curve.plot(256, ax=ax1, color=BLUE)
    color = ax1.lines[-1].get_color()
    add_patch(ax1, curve._nodes, color)
    reduced.plot(256, ax=ax2, color=BLUE)
    color = ax2.lines[-1].get_color()
    add_patch(ax2, reduced._nodes, color)
    ax1.axis("scaled")
    ax2.axis("scaled")
    _plot_helpers.add_plot_boundary(ax2)
    save_image(figure, "curve_reduce.png")


def curve_reduce_approx(curve, reduced):
    """Image for :meth:`.curve.Curve.reduce` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    color = ax.lines[-1].get_color()
    add_patch(ax, curve._nodes, color, alpha=0.25, node_color=color)
    reduced.plot(256, ax=ax, color=GREEN)
    color = ax.lines[-1].get_color()
    add_patch(ax, reduced._nodes, color, alpha=0.25, node_color=color)
    ax.axis("scaled")
    _plot_helpers.add_plot_boundary(ax)
    save_image(ax.figure, "curve_reduce_approx.png")
