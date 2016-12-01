# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Helper to make images that are intended for docs."""


import os

from matplotlib import patches
from matplotlib import path as _path_mod
import matplotlib.pyplot as plt
import numpy as np
try:
    import seaborn  # pylint: disable=unused-import
except ImportError:
    pass

import bezier
from bezier import _intersection_helpers


_DOCS_DIR = os.path.abspath(os.path.dirname(__file__))
IMAGES_DIR = os.path.join(_DOCS_DIR, 'images')
NO_IMAGES = 'NO_IMAGES' in os.environ


def save_image(figure, filename):
    """Save an image to the docs images directory.

    Args:
        filename (str): The name of the file (not containing
            directory info).
    """
    path = os.path.join(IMAGES_DIR, filename)
    figure.savefig(path, bbox_inches='tight')
    plt.close(figure)


def linearization_error():
    """Image for :func:`.linearization_error` docstring."""
    curve = bezier.Curve(np.array([
        [0.0, 0.0],
        [3.0, 1.0],
        [9.0, -2.0],
    ]))
    line = bezier.Curve(curve._nodes[(0, -1), :])

    midpoints = np.vstack([
        curve.evaluate(0.5),
        line.evaluate(0.5),
    ])

    ax = curve.plot(256)
    line.plot(256, ax=ax)
    ax.plot(midpoints[:, 0], midpoints[:, 1],
            color='black', linestyle='dashed')

    ax.axis('scaled')
    save_image(ax.figure, 'linearization_error.png')


def newton_refine1():
    """Image for :func:`.newton_refine` docstring."""
    curve1 = bezier.Curve(np.array([
        [0.0, 0.0],
        [2.0, 4.0],
        [4.0, 0.0],
    ]))
    curve2 = bezier.Curve(np.array([
        [2.0, 0.0],
        [0.0, 3.0],
    ]))

    # We are "interested" in the two incorrect s and t values
    # and the values they move to.
    s = 0.375
    new_s = s - 9.0 / 64.0
    t = 0.25
    new_t = t + 9.0 / 32.0
    points = np.vstack([
        curve1.evaluate(s),
        curve2.evaluate(t),
    ])
    points_new = np.vstack([
        curve1.evaluate(new_s),
        curve2.evaluate(new_t),
    ])

    ax = curve1.plot(256)
    curve2.plot(256, ax=ax)
    ax.plot(points[:, 0], points[:, 1],
            color='black', linestyle='None', marker='o',
            markeredgewidth=1, markerfacecolor='None')
    ax.plot(points_new[:, 0], points_new[:, 1],
            color='black', linestyle='None', marker='o')

    ax.axis('scaled')
    save_image(ax.figure, 'newton_refine1.png')


def newton_refine2():
    """Image for :func:`.newton_refine` docstring."""
    curve1 = bezier.Curve(np.array([
        [0.0, 0.0],
        [0.25, 2.0],
        [0.5, -2.0],
        [0.75, 2.0],
        [1.0, 0.0],
    ]))
    curve2 = bezier.Curve(np.array([
        [0.0, 1.0],
        [0.25, 0.5],
        [0.5, 0.5],
        [0.75, 0.5],
        [1.0, 0.0],
    ]))

    ax = curve1.plot(256)
    ax.lines[-1].zorder = 1
    curve2.plot(256, ax=ax)
    ax.lines[-1].zorder = 1

    s, t = 0.625, 0.625
    all_s = np.array([s, 0.0, 0.0, 0.0, 0.0])
    for index in range(4):
        s, t = _intersection_helpers.newton_refine(s, curve1, t, curve2)
        all_s[index + 1] = s

    points = curve1.evaluate_multi(all_s)
    colors = seaborn.dark_palette('blue', 5)
    ax.scatter(points[:, 0], points[:, 1], c=colors,
               alpha=0.75, zorder=2)

    ax.axis('scaled')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    save_image(ax.figure, 'newton_refine2.png')


def newton_refine3():
    """Image for :func:`.newton_refine` docstring."""
    curve1 = bezier.Curve(np.array([
        [0.0, 0.0],
        [0.5, 1.0],
        [1.0, 0.0],
    ]))
    curve2 = bezier.Curve(np.array([
        [0.0, 0.5],
        [1.0, 0.5],
    ]))

    ax = curve1.plot(256)
    ax.lines[-1].zorder = 1
    curve2.plot(256, ax=ax)
    ax.lines[-1].zorder = 1

    s, t = 0.375, 0.375
    all_s = np.array([s, 0.0, 0.0, 0.0, 0.0, 0.0])
    for index in range(5):
        s, t = _intersection_helpers.newton_refine(s, curve1, t, curve2)
        all_s[index + 1] = s

    points = curve1.evaluate_multi(all_s)
    colors = seaborn.dark_palette('blue', 6)
    ax.scatter(points[:, 0], points[:, 1], c=colors,
               s=20, alpha=0.75, zorder=2)

    ax.axis('scaled')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.5625)
    save_image(ax.figure, 'newton_refine3.png')


def segment_intersection1():
    """Image for :func:`.segment_intersection` docstring."""
    line0 = bezier.Curve(np.array([
        [0.0, 0.0],
        [2.0, 2.0],
    ]))
    line1 = bezier.Curve(np.array([
        [-1.0, 2.0],
        [1.0, 0.0],
    ]))

    ax = line0.plot(2)
    line1.plot(256, ax=ax)

    x_val, y_val = line0.evaluate(0.25)
    ax.plot([x_val], [y_val], color='black', marker='o')

    ax.axis('scaled')
    save_image(ax.figure, 'segment_intersection1.png')


def segment_intersection2():
    """Image for :func:`.segment_intersection` docstring."""
    line0 = bezier.Curve(np.array([
        [1.0, 0.0],
        [0.0, 1.0],
    ]))
    line1 = bezier.Curve(np.array([
        [-1.0, 3.0],
        [3.0, -1.0],
    ]))

    ax = line0.plot(2)
    line1.plot(2, ax=ax)

    ax.axis('scaled')
    save_image(ax.figure, 'segment_intersection2.png')


def helper_parallel_different(start0, end0, start1, end1, filename):
    figure = plt.figure()
    ax = figure.gca()

    points = np.vstack([start0, end0, start1, end1])
    ax.plot(points[:2, 0], points[:2, 1], marker='o')
    ax.plot(points[2:, 0], points[2:, 1], marker='o')

    ax.axis('scaled')
    left, bottom = np.min(points, axis=0)
    right, top = np.max(points, axis=0)
    width_x = right - left
    center_x = 0.5 * (right + left)
    width_y = top - bottom
    center_y = 0.5 * (top + bottom)

    ax.set_xlim(center_x - 0.5625 * width_x, center_x + 0.5625 * width_x)
    ax.set_ylim(center_y - 0.5625 * width_y, center_y + 0.5625 * width_y)

    save_image(figure, filename)


def parallel_different1():
    """Image for :func:`.parallel_different` docstring."""
    helper_parallel_different(
        np.array([[0.0, 1.0]]), np.array([[1.0, 1.0]]),
        np.array([[-1.0, 2.0]]), np.array([[3.0, 2.0]]),
        'parallel_different1.png')


def parallel_different2():
    """Image for :func:`.parallel_different` docstring."""
    helper_parallel_different(
        np.array([[1.0, 0.0]]), np.array([[3.0, -1.0]]),
        np.array([[4.0, -1.5]]), np.array([[5.0, -2.0]]),
        'parallel_different2.png')


def parallel_different3():
    """Image for :func:`.parallel_different` docstring."""
    helper_parallel_different(
        np.array([[1.0, 0.0]]), np.array([[3.0, -1.0]]),
        np.array([[-1.0, 1.0]]), np.array([[2.0, -0.5]]),
        'parallel_different3.png')


def parallel_different4():
    """Image for :func:`.parallel_different` docstring."""
    helper_parallel_different(
        np.array([[1.0, 0.0]]), np.array([[3.0, -1.0]]),
        np.array([[7.0, -3.0]]), np.array([[-3.0, 2.0]]),
        'parallel_different4.png')


def add_patch(ax, nodes, color, with_nodes=True):
    path = _path_mod.Path(nodes)
    patch = patches.PathPatch(
        path, facecolor=color, alpha=0.6)
    ax.add_patch(patch)
    if with_nodes:
        ax.plot(nodes[:, 0], nodes[:, 1], color='black',
                linestyle='None', marker='o')


def curve_constructor(curve):
    """Image for :class`.Curve` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256)
    line = ax.lines[0]

    nodes = curve._nodes
    ax.plot(nodes[:, 0], nodes[:, 1], color='black',
            linestyle='None', marker='o')
    add_patch(ax, nodes, line.get_color())

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'curve_constructor.png')


def curve_evaluate(curve):
    """Image for :meth`.Curve.evaluate` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256)
    points = curve.evaluate_multi(np.array([0.75]))
    ax.plot(points[:, 0], points[:, 1], color='black',
            linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'curve_evaluate.png')


def curve_subdivide(curve, left, right):
    """Image for :meth`.Curve.subdivide` docstring."""
    if NO_IMAGES:
        return

    figure = plt.figure()
    ax = figure.gca()
    add_patch(ax, curve._nodes, 'gray')

    ax = left.plot(256, ax=ax)
    line = ax.lines[-1]
    add_patch(ax, left._nodes, line.get_color())

    right.plot(256, ax=ax)
    line = ax.lines[-1]
    add_patch(ax, right._nodes, line.get_color())

    ax.axis('scaled')
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 3.125)
    save_image(ax.figure, 'curve_subdivide.png')


def curve_intersect(curve1, curve2, intersections):
    """Image for :meth`.Curve.intersect` docstring."""
    if NO_IMAGES:
        return

    ax = curve1.plot(256)
    curve2.plot(256, ax=ax)
    ax.plot(intersections[:, 0], intersections[:, 1],
            color='black', linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(0.0, 0.75)
    ax.set_ylim(0.0, 0.75)
    save_image(ax.figure, 'curve_intersect.png')


def surface_constructor():
    """Image for :class`.Surface` docstring."""
    nodes = np.array([
        [0.0, 0.0],
        [0.5, 0.0],
        [1.0, 0.25],
        [0.125, 0.5],
        [0.375, 0.375],
        [0.25, 1.0],
    ])
    surface = bezier.Surface(nodes)

    ax = surface.plot(256, with_nodes=True)
    line = ax.lines[0]
    add_patch(ax, nodes[(0, 1, 2, 5), :], line.get_color())
    delta = 1.0 / 32.0
    ax.text(nodes[0, 0], nodes[0, 1], r'$v_0$', fontsize=20,
            verticalalignment='top', horizontalalignment='right')
    ax.text(nodes[1, 0], nodes[1, 1], r'$v_1$', fontsize=20,
            verticalalignment='top', horizontalalignment='center')
    ax.text(nodes[2, 0], nodes[2, 1], r'$v_2$', fontsize=20,
            verticalalignment='top', horizontalalignment='left')
    ax.text(nodes[3, 0] - delta, nodes[3, 1], r'$v_3$', fontsize=20,
            verticalalignment='center', horizontalalignment='right')
    ax.text(nodes[4, 0] + delta, nodes[4, 1], r'$v_4$', fontsize=20,
            verticalalignment='center', horizontalalignment='left')
    ax.text(nodes[5, 0], nodes[5, 1] + delta, r'$v_5$', fontsize=20,
            verticalalignment='bottom', horizontalalignment='center')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, 'surface_constructor.png')


def surface_evaluate_barycentric():
    """Image for :meth`.Surface.evaluate_barycentric` docstring."""
    nodes = np.array([
        [0.0, 0.0],
        [0.5, 0.0],
        [1.0, 0.25],
        [0.125, 0.5],
        [0.375, 0.375],
        [0.25, 1.0],
    ])
    surface = bezier.Surface(nodes)

    ax = surface.plot(256)
    point = surface.evaluate_barycentric(0.125, 0.125, 0.75)
    ax.plot([point[0]], [point[1]], color='black',
            linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, 'surface_evaluate_barycentric.png')


def surface_evaluate_multi1():
    """Image for :meth`.Surface.evaluate_multi` docstring."""
    nodes = np.array([
        [0.0, 0.0],
        [2.0, 1.0],
        [-3.0, 2.0],
    ])
    surface = bezier.Surface(nodes)

    ax = surface.plot(256)
    param_vals = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 0.5],
    ])
    points = surface.evaluate_multi(param_vals)
    ax.plot(points[:, 0], points[:, 1], color='black',
            linestyle='None', marker='o')

    delta = 1.0 / 32.0
    ax.text(points[0, 0], points[0, 1], r'$w_0$', fontsize=20,
            verticalalignment='top', horizontalalignment='right')
    ax.text(points[1, 0] + 2 * delta, points[1, 1], r'$w_1$', fontsize=20,
            verticalalignment='center', horizontalalignment='left')
    ax.text(points[2, 0], points[2, 1] + delta, r'$w_2$', fontsize=20,
            verticalalignment='bottom', horizontalalignment='left')

    ax.axis('scaled')
    ax.set_xlim(-3.125, 2.375)
    ax.set_ylim(-0.25, 2.125)
    save_image(ax.figure, 'surface_evaluate_multi1.png')


def surface_evaluate_multi2():
    """Image for :meth`.Surface.evaluate_multi` docstring."""
    surface = bezier.Surface(np.array([
        [0.0, 0.0],
        [1.0, 0.75],
        [2.0, 1.0],
        [-1.5, 1.0],
        [-0.5, 1.5],
        [-3.0, 2.0],
    ]))

    ax = surface.plot(256)
    param_vals = np.array([
        [0.0, 0.25, 0.75],
        [1.0, 0.0, 0.0],
        [0.25, 0.5, 0.25],
        [0.375, 0.25, 0.375],
    ])
    points = surface.evaluate_multi(param_vals)
    ax.plot(points[:, 0], points[:, 1], color='black',
            linestyle='None', marker='o')

    delta = 1.0 / 32.0
    ax.text(points[0, 0], points[0, 1] + delta, r'$w_0$', fontsize=20,
            verticalalignment='bottom', horizontalalignment='center')
    ax.text(points[1, 0], points[1, 1] - delta, r'$w_1$', fontsize=20,
            verticalalignment='top', horizontalalignment='right')
    ax.text(points[2, 0], points[2, 1], r'$w_2$', fontsize=20,
            verticalalignment='bottom', horizontalalignment='left')
    ax.text(points[3, 0], points[3, 1], r'$w_3$', fontsize=20,
            verticalalignment='top', horizontalalignment='right')

    ax.axis('scaled')
    ax.set_xlim(-3.125, 2.125)
    ax.set_ylim(-0.3125, 2.125)
    save_image(ax.figure, 'surface_evaluate_multi2.png')


def surface_is_valid1():
    """Image for :meth`.Surface.is_valid` docstring."""
    surface = bezier.Surface(np.array([
        [0.0, 0.0],
        [1.0, 1.0],
        [2.0, 2.0],
    ]))

    ax = surface.plot(256)
    ax.axis('scaled')
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 2.125)
    save_image(ax.figure, 'surface_is_valid1.png')


def surface_is_valid2():
    """Image for :meth`.Surface.is_valid` docstring."""
    surface = bezier.Surface(np.array([
        [0.0, 0.0],
        [0.5, 0.125],
        [1.0, 0.0],
        [-0.125, 0.5],
        [0.5, 0.5],
        [0.0, 1.0],
    ]))

    ax = surface.plot(256)
    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.0625)
    ax.set_ylim(-0.0625, 1.0625)
    save_image(ax.figure, 'surface_is_valid2.png')


def surface_is_valid3():
    """Image for :meth`.Surface.is_valid` docstring."""
    surface = bezier.Surface(np.array([
        [1.0, 0.0],
        [0.0, 0.0],
        [1.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 1.0],
    ]))
    edge1, edge2, edge3 = surface.edges

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
    jacobian_zero_params = np.zeros((2 * N - 1, 2))
    jacobian_zero_params[:N, 0] = s_vals
    jacobian_zero_params[:N, 1] = t_top
    jacobian_zero_params[N:, 0] = s_vals[-2::-1]
    jacobian_zero_params[N:, 1] = t_bottom[-2::-1]
    jac_edge = surface.evaluate_multi(jacobian_zero_params)

    # Add the surface to the plot and add a dashed line
    # for each "true" edge.
    figure = plt.figure()
    ax = figure.gca()
    line, = ax.plot(jac_edge[:, 0], jac_edge[:, 1])
    color = line.get_color()

    ax.plot(points1[:, 0], points1[:, 1],
            color='black', linestyle='dashed')
    ax.plot(points2[:, 0], points2[:, 1],
            color='black', linestyle='dashed')
    ax.plot(points3[:, 0], points3[:, 1],
            color='black', linestyle='dashed')

    polygon = np.vstack([
        points1[1:, :],
        points2[1:, :],
        jac_edge[1:, :],
    ])
    add_patch(ax, polygon, color, with_nodes=False)

    ax.axis('scaled')
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.0625, 1.0625)
    save_image(ax.figure, 'surface_is_valid3.png')


def surface_subdivide1():
    """Image for :meth`.Surface.subdivide` docstring."""
    surface = bezier.Surface(np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ]))
    surf_a, surf_b, surf_c, surf_d = surface.subdivide()

    figure, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    for ax in (ax1, ax2, ax3, ax4):
        surface.plot(2, ax=ax)

    surf_a.plot(2, ax=ax1)
    ax1.text(1.0 / 6.0, 1.0 / 6.0, r'$A$', fontsize=20,
             verticalalignment='center', horizontalalignment='center')

    surf_b.plot(2, ax=ax2)
    ax2.text(1.0 / 3.0, 1.0 / 3.0, r'$B$', fontsize=20,
             verticalalignment='center', horizontalalignment='center')

    surf_c.plot(2, ax=ax3)
    ax3.text(2.0 / 3.0, 1.0 / 6.0, r'$C$', fontsize=20,
             verticalalignment='center', horizontalalignment='center')

    surf_d.plot(2, ax=ax4)
    ax4.text(1.0 / 6.0, 2.0 / 3.0, r'$D$', fontsize=20,
             verticalalignment='center', horizontalalignment='center')

    for ax in (ax1, ax2, ax3, ax4):
        ax.axis('scaled')

    save_image(figure, 'surface_subdivide1')


def add_edges(ax, surface, s_vals, color):
    edge1, edge2, edge3 = surface.edges
    # Compute points on each edge.
    points1 = edge1.evaluate_multi(s_vals)
    points2 = edge2.evaluate_multi(s_vals)
    points3 = edge3.evaluate_multi(s_vals)

    # Add the points to the plot.
    ax.plot(points1[:, 0], points1[:, 1], color=color)
    ax.plot(points2[:, 0], points2[:, 1], color=color)
    ax.plot(points3[:, 0], points3[:, 1], color=color)


def surface_subdivide2():
    """Image for :meth`.Surface.subdivide` docstring."""
    # Plot set-up.
    figure = plt.figure()
    ax = figure.gca()
    colors = seaborn.husl_palette(6)

    # Define surface and sub-surface.
    surface = bezier.Surface(np.array([
        [-1.0, 0.0],
        [0.5, 0.5],
        [2.0, 0.0],
        [0.25, 1.75],
        [2.0, 3.0],
        [0.0, 4.0],
    ]))
    _, surf_b, _, _ = surface.subdivide()

    N = 128
    s_vals = np.linspace(0.0, 1.0, N + 1)
    # Add edges from surface.
    add_edges(ax, surface, s_vals, colors[4])
    # Now do the same for surface B.
    add_edges(ax, surf_b, s_vals, colors[0])

    # Add the control points polygon for the original surface.
    nodes = surface._nodes[(0, 2, 4, 5, 0), :]
    add_patch(ax, nodes, colors[2], with_nodes=False)

    # Add the control points polygon for the sub-surface.
    nodes = surf_b._nodes[(0, 1, 2, 5, 3, 0), :]
    add_patch(ax, nodes, colors[1])
    # Take those same points and add the boundary.
    ax.plot(nodes[:, 0], nodes[:, 1],
            color='black', linestyle='dashed')

    ax.axis('scaled')
    ax.set_xlim(-1.125, 2.125)
    ax.set_ylim(-0.125, 4.125)
    save_image(ax.figure, 'surface_subdivide2')


def curved_polygon_constructor1(curved_poly):
    """Image for :class`.CurvedPolygon` docstring."""
    if NO_IMAGES:
        return

    ax = curved_poly.plot(256)
    ax.axis('scaled')
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.625, 1.625)
    save_image(ax.figure, 'curved_polygon_constructor1.png')


def curved_polygon_constructor2(curved_poly):
    """Image for :class`.CurvedPolygon` docstring."""
    if NO_IMAGES:
        return

    ax = curved_poly.plot(256)
    ax.axis('scaled')
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, 'curved_polygon_constructor2.png')


def main():
    linearization_error()
    newton_refine1()
    newton_refine2()
    newton_refine3()
    segment_intersection1()
    segment_intersection2()
    parallel_different1()
    parallel_different2()
    parallel_different3()
    parallel_different4()
    surface_constructor()
    surface_evaluate_barycentric()
    surface_evaluate_multi1()
    surface_evaluate_multi2()
    surface_is_valid1()
    surface_is_valid2()
    surface_is_valid3()
    surface_subdivide1()
    surface_subdivide2()


if __name__ == '__main__':
    main()
