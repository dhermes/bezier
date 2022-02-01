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

   $ nox --session docs_images
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
from bezier import _geometric_intersection
from bezier import _helpers
from bezier import _plot_helpers
from bezier.hazmat import clipping
from bezier.hazmat import geometric_intersection as _py_geometric_intersection


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


def curve_self_intersect2(curve, self_intersections):
    """Image for :meth`.Curve.self_intersections` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    if self_intersections.shape != (2, 1):
        raise ValueError("Unexpected shape", self_intersections)

    s1_val = self_intersections[0, 0]
    intersection_xy = curve.evaluate(s1_val)
    ax.plot(
        intersection_xy[0, :],
        intersection_xy[1, :],
        color="black",
        linestyle="None",
        marker="o",
    )
    ax.axis("scaled")
    ax.set_xlim(-0.8125, 0.0625)
    ax.set_ylim(0.75, 2.125)
    save_image(ax.figure, "curve_self_intersect2.png")


def curve_self_intersect3(curve, self_intersections):
    """Image for :meth`.Curve.self_intersections` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256, color=BLUE)
    if self_intersections.shape != (2, 2):
        raise ValueError("Unexpected shape", self_intersections)

    s1_vals = np.asfortranarray(self_intersections[0, :])
    intersection_xy = curve.evaluate_multi(s1_vals)
    ax.plot(
        intersection_xy[0, :],
        intersection_xy[1, :],
        color="black",
        linestyle="None",
        marker="o",
    )
    ax.axis("scaled")
    ax.set_xlim(-330.0, 330.0)
    ax.set_ylim(0.125, 266.0)
    save_image(ax.figure, "curve_self_intersect3.png")


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


def simple_axis(ax):
    ax.axis("scaled")
    ax.set_xticklabels([])
    ax.set_yticklabels([])


def plot_with_bbox(curve, ax, color, with_nodes=False):
    curve.plot(256, color=color, ax=ax)
    left, right, bottom, top = _helpers.bbox(curve._nodes)
    bbox_nodes = np.asfortranarray(
        [[left, right, right, left], [bottom, bottom, top, top]]
    )
    add_patch(ax, bbox_nodes, color, with_nodes=False)
    if with_nodes:
        ax.plot(
            curve._nodes[0, :],
            curve._nodes[1, :],
            color=color,
            linestyle="None",
            marker="o",
            markersize=4,
        )


def plot_with_convex_hull(curve, ax, color, with_nodes=False):
    curve.plot(256, color=color, ax=ax)
    convex_hull = _helpers.simple_convex_hull(curve._nodes)
    add_patch(ax, convex_hull, color, with_nodes=False)
    if with_nodes:
        ax.plot(
            curve._nodes[0, :],
            curve._nodes[1, :],
            color=color,
            linestyle="None",
            marker="o",
            markersize=4,
        )


def _curve_boundary_predicate(filename, curve_boundary_plot, headers):
    header1, header2, header3 = headers

    figure, (ax1, ax2, ax3) = plt.subplots(1, 3)

    control_pts1a = np.asfortranarray([[0.0, 0.375, 1.0], [0.0, 0.5, 0.125]])
    curve1a = bezier.Curve(control_pts1a, degree=2)
    control_pts1b = np.asfortranarray(
        [[0.25, -0.125, 0.5], [-0.125, 0.375, 1.0]]
    )
    curve1b = bezier.Curve(control_pts1b, degree=2)
    curve_boundary_plot(curve1a, ax1, BLUE)
    curve_boundary_plot(curve1b, ax1, GREEN)

    control_pts2a = np.asfortranarray([[0.0, 0.75, 1.0], [1.0, 0.75, 0.0]])
    curve2a = bezier.Curve(control_pts2a, degree=2)
    control_pts2b = np.asfortranarray(
        [[0.625, 0.875, 1.625], [1.625, 0.875, 0.625]]
    )
    curve2b = bezier.Curve(control_pts2b, degree=2)
    curve_boundary_plot(curve2a, ax2, BLUE)
    curve_boundary_plot(curve2b, ax2, GREEN)

    control_pts3a = np.asfortranarray([[0.0, 0.25, 1.0], [-0.25, 0.25, -0.75]])
    curve3a = bezier.Curve(control_pts3a, degree=2)
    control_pts3b = np.asfortranarray([[1.0, 1.5, 2.0], [-1.0, -1.5, -1.0]])
    curve3b = bezier.Curve(control_pts3b, degree=2)
    curve_boundary_plot(curve3a, ax3, BLUE)
    curve_boundary_plot(curve3b, ax3, GREEN)

    for ax in (ax1, ax2, ax3):
        simple_axis(ax)

    text_size = 10
    ax1.set_xlim(-0.2, 1.1)
    ax1.set_ylim(-0.2, 1.1)
    ax1.set_title(header1, fontsize=text_size)
    ax2.set_xlim(-0.1, 1.75)
    ax2.set_ylim(-0.1, 1.75)
    ax2.set_title(header2, fontsize=text_size)
    ax3.set_xlim(-0.1, 2.1)
    ax3.set_ylim(-1.7, 0.5)
    ax3.set_title(header3, fontsize=text_size)

    figure.set_size_inches(6.0, 2.2)
    figure.subplots_adjust(
        left=0.01, bottom=0.01, right=0.99, top=0.9, wspace=0.04, hspace=0.2
    )
    save_image(figure, filename)


def bounding_box_predicate():
    headers = ("MAYBE", "MAYBE", "NO")
    _curve_boundary_predicate(
        "bounding_box_predicate.png", plot_with_bbox, headers
    )


def convex_hull_predicate():
    headers = ("MAYBE", "NO", "NO")
    _curve_boundary_predicate(
        "convex_hull_predicate.png", plot_with_convex_hull, headers
    )


def subdivide_curve():
    figure, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True)
    nodes = np.asfortranarray([[0.0, 1.0, 2.0, 4.0], [0.0, 4.0, 0.0, 3.0]])
    curve = bezier.Curve.from_nodes(nodes)
    left, right = curve.subdivide()
    curve.plot(256, ax=ax1, alpha=0.25, color="black")
    left.plot(256, ax=ax1)
    curve.plot(256, ax=ax2)
    curve.plot(256, ax=ax3, alpha=0.25, color="black")
    right.plot(256, ax=ax3)
    text_size = 10
    ax1.text(
        2.5,
        0.25,
        r"$\left[0, \frac{1}{2}\right]$",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=text_size,
    )
    ax2.text(
        2.5,
        0.25,
        r"$\left[0, 1\right]$",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=text_size,
    )
    ax3.text(
        2.5,
        0.25,
        r"$\left[\frac{1}{2}, 1\right]$",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=text_size,
    )

    for ax in (ax1, ax2, ax3):
        simple_axis(ax)

    figure.set_size_inches(6.0, 1.5)
    figure.subplots_adjust(
        left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0.04, hspace=0.2
    )
    save_image(figure, "subdivide_curve.png")


def bbox_intersect(curve1, curve2):
    enum_val = _geometric_intersection.bbox_intersect(
        curve1.nodes, curve2.nodes
    )
    return enum_val != _py_geometric_intersection.BoxIntersectionType.DISJOINT


def refine_candidates(left, right):
    new_left = []
    for curve in left:
        new_left.extend(curve.subdivide())

    new_right = []
    for curve in right:
        new_right.extend(curve.subdivide())

    keep_left = []
    keep_right = []
    for curve1 in new_left:
        for curve2 in new_right:
            if bbox_intersect(curve1, curve2):
                keep_left.append(curve1)
                if curve2 not in keep_right:
                    keep_right.append(curve2)

    return keep_left, keep_right


def unique_curves(pairs):
    left_tuples = set()
    right_tuples = set()

    left_curves = []
    right_curves = []
    for left, right in pairs:
        as_tuple = tuple(left._nodes.flatten(order="F"))
        if as_tuple not in left_tuples:
            left_tuples.add(as_tuple)
            left_curves.append(left)

        as_tuple = tuple(right._nodes.flatten(order="F"))
        if as_tuple not in right_tuples:
            right_tuples.add(as_tuple)
            right_curves.append(right)

    return left_curves, right_curves


def subdivision_process():
    nodes15 = np.asfortranarray([[0.25, 0.625, 1.0], [0.625, 0.25, 1.0]])
    curve15 = bezier.Curve(nodes15, degree=2)
    nodes25 = np.asfortranarray([[0.0, 0.25, 0.75, 1.0], [0.5, 1.0, 1.5, 0.5]])
    curve25 = bezier.Curve(nodes25, degree=3)

    figure, all_axes = plt.subplots(2, 3, sharex=True, sharey=True)
    ax1, ax2, ax3, ax4, ax5, ax6 = all_axes.flatten()

    color1 = BLUE
    color2 = GREEN
    plot_with_bbox(curve15, ax1, color1)
    plot_with_bbox(curve25, ax1, color2)

    left, right = refine_candidates([curve15], [curve25])
    for curve in left:
        plot_with_bbox(curve, ax2, color1)
    for curve in right:
        plot_with_bbox(curve, ax2, color2)

    for ax in (ax3, ax4, ax5, ax6):
        left, right = refine_candidates(left, right)
        curve15.plot(256, color=color1, alpha=0.5, ax=ax)
        for curve in left:
            plot_with_bbox(curve, ax, color=color1)
        curve25.plot(256, color=color2, alpha=0.5, ax=ax)
        for curve in right:
            plot_with_bbox(curve, ax, color2)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
        simple_axis(ax)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0.4, 1.15)

    figure.set_size_inches(6.0, 2.8)
    figure.subplots_adjust(
        left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0.04, hspace=0.04
    )
    save_image(figure, "subdivision_process.png")


def _subdivision_pruning_zoom(all_axes, column):
    half_width = 0.5 ** (column + 2)
    min_x = 0.75 - half_width
    max_x = 0.75 + half_width
    min_y = 0.25 - half_width
    max_y = 0.25 + half_width

    for row in (0, 1, 2):
        all_axes[row, column].plot(
            [min_x, max_x, max_x, min_x, min_x],
            [min_y, min_y, max_y, max_y, min_y],
            color="black",
        )

    buffer = 0.5 ** (column + 6)
    for row in (1, 2):
        all_axes[row, column].set_xlim(min_x - buffer, max_x + buffer)
        all_axes[row, column].set_ylim(min_y - buffer, max_y + buffer)


def subdivision_pruning():
    figure, all_axes = plt.subplots(3, 4)

    nodes69 = np.asfortranarray([[0.0, 1.0, 1.0], [0.0, 0.0, 1.0]])
    curve69 = bezier.Curve(nodes69, degree=2)
    delta = np.asfortranarray([[1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]]) / 32.0
    nodes_other = nodes69 + delta
    curve_other = bezier.Curve(nodes_other, degree=2)

    color1 = BLUE
    color2 = GREEN
    for ax in all_axes.flatten():
        curve69.plot(256, color=color1, ax=ax)
        curve_other.plot(256, color=color2, ax=ax)

    candidates = {0: [(curve69, curve_other)]}
    intersections = []
    for i in range(5):
        candidates[i + 1] = _py_geometric_intersection.intersect_one_round(
            candidates[i], intersections
        )

    for column in (0, 1, 2, 3):
        left_curves, right_curves = unique_curves(candidates[column + 2])
        for curve in left_curves:
            plot_with_bbox(curve, all_axes[0, column], color1)
            plot_with_bbox(curve, all_axes[1, column], color1, with_nodes=True)
            plot_with_convex_hull(
                curve, all_axes[2, column], color1, with_nodes=True
            )
        for curve in right_curves:
            plot_with_bbox(curve, all_axes[0, column], color2)
            plot_with_bbox(curve, all_axes[1, column], color2, with_nodes=True)
            plot_with_convex_hull(
                curve, all_axes[2, column], color2, with_nodes=True
            )

    for ax in all_axes.flatten():
        simple_axis(ax)

    _subdivision_pruning_zoom(all_axes, 0)
    _subdivision_pruning_zoom(all_axes, 1)
    _subdivision_pruning_zoom(all_axes, 2)
    _subdivision_pruning_zoom(all_axes, 3)

    intersection_params = curve69.intersect(curve_other)
    s_vals = intersection_params[0, :]
    intersections = curve69.evaluate_multi(s_vals)
    for column in (0, 1, 2, 3):
        all_axes[0, column].plot(
            intersections[0, :],
            intersections[1, :],
            color="black",
            linestyle="None",
            marker="o",
            markersize=4,
        )

    save_image(figure, "subdivision_pruning.png")


def _plot_endpoints_line(ax, fat_line_coeffs, **plot_kwargs):
    # NOTE: This assumes the x-limits have already been set for the axis.
    coeff_a, coeff_b, coeff_c, _, _ = fat_line_coeffs
    if coeff_b == 0.0:
        raise NotImplementedError("Vertical lines not supported")

    min_x, max_x = ax.get_xlim()
    # ax + by + c = 0 ==> y = -(ax + c)/b
    min_y = -(coeff_a * min_x + coeff_c) / coeff_b
    max_y = -(coeff_a * max_x + coeff_c) / coeff_b

    ax.plot(
        [min_x, max_x],
        [min_y, max_y],
        **plot_kwargs,
    )
    ax.set_xlim(min_x, max_x)


def _normalize_implicit_line_tuple(info):
    length = np.linalg.norm(info[:2], ord=2)
    return tuple(np.array(info) / length)


def compute_implicit_line(nodes):
    """Image for :func:`.hazmat.clipping.compute_implicit_line` docstring."""
    if NO_IMAGES:
        return

    curve = bezier.Curve.from_nodes(nodes)
    ax = curve.plot(256, color=BLUE)

    min_x, max_x = nodes[0, (0, -1)]
    min_y, max_y = nodes[1, (0, -1)]
    ax.plot(
        [min_x, max_x, max_x],
        [min_y, min_y, max_y],
        color="black",
        linestyle="dashed",
    )

    ax.axis("scaled")
    # NOTE: This "cheats" and assumes knowledge of what's actually in ``nodes``.
    ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0])
    ax.set_yticks([0.0, 1.0, 2.0, 3.0])

    info = clipping.compute_fat_line(nodes)
    info = _normalize_implicit_line_tuple(info)
    _plot_endpoints_line(ax, info, color="black")
    save_image(ax.figure, "compute_implicit_line.png")


def _plot_fat_lines(ax, fat_line_coeffs, **fill_between_kwargs):
    # NOTE: This assumes the x-limits have already been set for the axis.
    coeff_a, coeff_b, coeff_c, d_low, d_high = fat_line_coeffs
    if coeff_b == 0.0:
        raise NotImplementedError("Vertical lines not supported")

    min_x, max_x = ax.get_xlim()
    coeff_c_low = coeff_c - d_low
    coeff_c_high = coeff_c - d_high
    # ax + by + c = 0 ==> y = -(ax + c)/b
    min_y_low = -(coeff_a * min_x + coeff_c_low) / coeff_b
    min_y_high = -(coeff_a * min_x + coeff_c_high) / coeff_b
    max_y_low = -(coeff_a * max_x + coeff_c_low) / coeff_b
    max_y_high = -(coeff_a * max_x + coeff_c_high) / coeff_b

    ax.fill_between(
        [min_x, max_x],
        [min_y_low, max_y_low],
        [min_y_high, max_y_high],
        **fill_between_kwargs,
    )
    ax.set_xlim(min_x, max_x)


def _add_perpendicular_segments(ax, nodes, fat_line_coeffs, color):
    coeff_a, coeff_b, coeff_c, _, _ = fat_line_coeffs

    _, num_nodes = nodes.shape
    for index in range(num_nodes):
        # ax + by + c = 0 is perpendicular to lines of the form
        # bx - ay = c'
        curr_x, curr_y = nodes[:, index]
        c_prime = coeff_b * curr_x - coeff_a * curr_y
        # bx - ay = c' intersects ax + by + c = 0 at
        # [x0, y0] = [b c' - a c, -a c' - b c] (assuming a^2 + b^2 == 1)
        x_intersect = coeff_b * c_prime - coeff_a * coeff_c
        y_intersect = -coeff_a * c_prime - coeff_b * coeff_c
        ax.plot(
            [curr_x, x_intersect],
            [curr_y, y_intersect],
            color=color,
            linestyle="dashed",
        )


def compute_fat_line(nodes, fat_line_coeffs):
    """Image for :func:`.hazmat.clipping.compute_fat_line` docstring."""
    if NO_IMAGES:
        return

    fat_line_coeffs = _normalize_implicit_line_tuple(fat_line_coeffs)
    curve = bezier.Curve.from_nodes(nodes)
    ax = curve.plot(256, color=BLUE)
    ax.plot(
        nodes[0, :],
        nodes[1, :],
        marker="o",
        color=BLUE,
        linestyle="none",
    )
    _add_perpendicular_segments(ax, nodes, fat_line_coeffs, BLUE)

    ax.axis("scaled")
    _plot_endpoints_line(ax, fat_line_coeffs, color=BLUE, linestyle="dashed")
    _plot_fat_lines(ax, fat_line_coeffs, color=BLUE, alpha=0.5)
    save_image(ax.figure, "compute_fat_line.png")


def clip_range(nodes1, nodes2):
    """Image for :func:`.hazmat.clipping.clip_range` docstring."""
    if NO_IMAGES:
        return

    curve1 = bezier.Curve.from_nodes(nodes1)
    curve2 = bezier.Curve.from_nodes(nodes2)

    # Plot both curves as well as the nodes.
    ax = curve1.plot(256, color=BLUE)
    curve2.plot(256, ax=ax, color=GREEN)
    ax.plot(
        nodes1[0, :],
        nodes1[1, :],
        marker="o",
        color=BLUE,
        linestyle="none",
    )
    ax.plot(
        nodes2[0, :],
        nodes2[1, :],
        marker="o",
        color=GREEN,
        linestyle="none",
    )

    fat_line_coeffs = clipping.compute_fat_line(nodes1)
    fat_line_coeffs = _normalize_implicit_line_tuple(fat_line_coeffs)
    # Add perpendicular lines to the "implicit" line.
    _add_perpendicular_segments(ax, nodes2, fat_line_coeffs, GREEN)

    # Establish boundary **assuming** contents of ``nodes1`` and ``nodes2``.
    ax.axis("scaled")
    ax.set_xlim(-0.625, 7.375)
    ax.set_ylim(-0.25, 4.625)

    _plot_endpoints_line(ax, fat_line_coeffs, color=BLUE, linestyle="dashed")
    _plot_fat_lines(ax, fat_line_coeffs, color=BLUE, alpha=0.5)

    save_image(ax.figure, "clip_range.png")


def clip_range_distances(nodes1, nodes2):
    """Image for :func:`.hazmat.clipping.clip_range` docstring."""
    if NO_IMAGES:
        return

    figure = plt.figure()
    ax = figure.gca()

    fat_line_coeffs = clipping.compute_fat_line(nodes1)
    coeff_a, coeff_b, coeff_c, d_min, d_max = fat_line_coeffs
    degree2, polynomial = clipping._clip_range_polynomial(
        nodes2, coeff_a, coeff_b, coeff_c
    )
    ax.fill_between([0.0, degree2], d_min, d_max, color=BLUE, alpha=0.25)
    s_min, s_max = clipping.clip_range(nodes1, nodes2)

    convex_hull = _helpers.simple_convex_hull(polynomial)
    add_patch(
        ax,
        convex_hull,
        GREEN,
        with_nodes=True,
        alpha=0.625,
        node_color=GREEN,
    )

    # Plot the true distance function ``d(t)``.
    t_values = np.linspace(0.0, 1.0, 257)
    curve2 = bezier.Curve.from_nodes(nodes2)
    evaluated = curve2.evaluate_multi(t_values)
    x_values = degree2 * t_values
    d_values = coeff_a * evaluated[0, :] + coeff_b * evaluated[1, :] + coeff_c
    ax.plot(x_values, d_values, color=GREEN)

    # Add dashed lines to each control point in the convex hull.
    for index in range(degree2 + 1):
        x_val, y_val = polynomial[:, index]
        ax.plot([x_val, x_val], [0.0, y_val], color=GREEN, linestyle="dashed")

    # NOTE: This "cheats" and uses the fact that it knows that ``s_min``
    #       corresponds to ``d_max`` and ``s_max`` corresponds to ``d_min``.
    ax.plot(
        [degree2 * s_min, degree2 * s_max],
        [d_max, d_min],
        color="black",
        marker="o",
        linestyle="none",
    )

    # Use minor xticks **above** for showing s_min and s_max.
    jitter = 0.5**5
    # NOTE: We introduce ``jitter`` to avoid using the same value for a minor
    #       xtick that is used for a major one. When ``matplotlib`` sees a
    #       minor xtick at the exact same value used by a major xtick, it
    #       ignores the tick.
    ax.set_xticks(
        [degree2 * s_min + jitter, degree2 * s_max - jitter], minor=True
    )
    ax.set_xticklabels([f"$t = {s_min}$", f"$t = {s_max}$"], minor=True)
    ax.tick_params(
        axis="x",
        which="minor",
        direction="in",
        top=False,
        bottom=False,
        labelbottom=False,
        labeltop=True,
    )
    # Add line up to minor xticks. Similar to the dots on ``s_min`` and
    # ``s_max`` this "cheats" with the correspondence to ``d_min`` / ``d_max``.
    min_y, max_y = ax.get_ylim()
    ax.plot(
        [degree2 * s_min, degree2 * s_min],
        [d_max, max_y],
        color="black",
        alpha=0.125,
        linestyle="dashed",
    )
    ax.plot(
        [degree2 * s_max, degree2 * s_max],
        [d_min, max_y],
        color="black",
        alpha=0.125,
        linestyle="dashed",
    )
    ax.set_ylim(min_y, max_y)

    ax.set_xlabel("$2t$")
    ax.set_ylabel("$d(t)$", rotation=0)
    save_image(figure, "clip_range_distances.png")


def main():
    bounding_box_predicate()
    convex_hull_predicate()
    subdivide_curve()
    subdivision_process()
    subdivision_pruning()


if __name__ == "__main__":
    main()
