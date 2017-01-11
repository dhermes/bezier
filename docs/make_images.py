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

"""Helper to make images that are intended for docs.

To actually execute these functions with the desired inputs, run:

.. code-block:: console

   $ tox -e docs-images
"""


import fractions
import os

from matplotlib import patches
from matplotlib import path as _path_mod
import matplotlib.pyplot as plt
import numpy as np
try:
    import seaborn
except ImportError:
    pass

import bezier
from bezier import _intersection_helpers
from bezier import _plot_helpers


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


def linearization_error(curve):
    """Image for :func:`.linearization_error` docstring."""
    if NO_IMAGES:
        return

    line = bezier.Curve.from_nodes(curve._nodes[(0, -1), :])
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


def newton_refine1(s, new_s, curve1, t, new_t, curve2):
    """Image for :func:`.newton_refine` docstring."""
    if NO_IMAGES:
        return

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


def newton_refine2(s_vals, curve1, curve2):
    """Image for :func:`.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve1.plot(256)
    ax.lines[-1].zorder = 1
    curve2.plot(256, ax=ax)
    ax.lines[-1].zorder = 1

    points = curve1.evaluate_multi(np.array(s_vals))
    colors = seaborn.dark_palette('blue', 5)
    ax.scatter(points[:, 0], points[:, 1], c=colors,
               s=20, alpha=0.75, zorder=2)

    ax.axis('scaled')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    save_image(ax.figure, 'newton_refine2.png')


def newton_refine3(s_vals, curve1, curve2):
    """Image for :func:`.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve1.plot(256)
    ax.lines[-1].zorder = 1
    curve2.plot(256, ax=ax)
    ax.lines[-1].zorder = 1

    points = curve1.evaluate_multi(np.array(s_vals))
    colors = seaborn.dark_palette('blue', 6)
    ax.scatter(points[:, 0], points[:, 1], c=colors,
               s=20, alpha=0.75, zorder=2)

    ax.axis('scaled')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.5625)
    save_image(ax.figure, 'newton_refine3.png')


def segment_intersection1(start0, end0, start1, end1, s):
    """Image for :func:`.segment_intersection` docstring."""
    if NO_IMAGES:
        return

    line0 = bezier.Curve.from_nodes(np.vstack([start0, end0]))
    line1 = bezier.Curve.from_nodes(np.vstack([start1, end1]))

    ax = line0.plot(2)
    line1.plot(256, ax=ax)

    (x_val, y_val), = line0.evaluate(s)
    ax.plot([x_val], [y_val], color='black', marker='o')

    ax.axis('scaled')
    save_image(ax.figure, 'segment_intersection1.png')


def segment_intersection2(start0, end0, start1, end1):
    """Image for :func:`.segment_intersection` docstring."""
    if NO_IMAGES:
        return

    line0 = bezier.Curve.from_nodes(np.vstack([start0, end0]))
    line1 = bezier.Curve.from_nodes(np.vstack([start1, end1]))

    ax = line0.plot(2)
    line1.plot(2, ax=ax)

    ax.axis('scaled')
    save_image(ax.figure, 'segment_intersection2.png')


def helper_parallel_different(start0, end0, start1, end1, filename):
    """Image for :func:`.parallel_different` docstring."""
    if NO_IMAGES:
        return

    figure = plt.figure()
    ax = figure.gca()

    points = np.vstack([start0, end0, start1, end1])
    ax.plot(points[:2, 0], points[:2, 1], marker='o')
    ax.plot(points[2:, 0], points[2:, 1], marker='o')

    ax.axis('scaled')
    _plot_helpers.add_plot_boundary(ax)

    save_image(figure, filename)


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


def surface_constructor(surface):
    """Image for :class`.Surface` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256, with_nodes=True)
    line = ax.lines[0]

    nodes = surface._nodes
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


def surface_evaluate_barycentric(surface, point):
    """Image for :meth`.Surface.evaluate_barycentric` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256)
    ax.plot(point[:, 0], point[:, 1], color='black',
            linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, 'surface_evaluate_barycentric.png')


def surface_evaluate_multi1(surface, points):
    """Image for :meth`.Surface.evaluate_multi` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256)
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


def surface_evaluate_multi2(surface, points):
    """Image for :meth`.Surface.evaluate_multi` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256)
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


def surface_is_valid1(surface):
    """Image for :meth`.Surface.is_valid` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256)
    ax.axis('scaled')
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 2.125)
    save_image(ax.figure, 'surface_is_valid1.png')


def surface_is_valid2(surface):
    """Image for :meth`.Surface.is_valid` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256)
    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.0625)
    ax.set_ylim(-0.0625, 1.0625)
    save_image(ax.figure, 'surface_is_valid2.png')


def surface_is_valid3(surface):
    """Image for :meth`.Surface.is_valid` docstring."""
    if NO_IMAGES:
        return

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
    if NO_IMAGES:
        return

    surface = bezier.Surface.from_nodes(np.array([
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


def surface_subdivide2(surface, sub_surface_b):
    """Image for :meth`.Surface.subdivide` docstring."""
    if NO_IMAGES:
        return

    # Plot set-up.
    figure = plt.figure()
    ax = figure.gca()
    colors = seaborn.husl_palette(6)

    N = 128
    s_vals = np.linspace(0.0, 1.0, N + 1)
    # Add edges from surface.
    add_edges(ax, surface, s_vals, colors[4])
    # Now do the same for surface B.
    add_edges(ax, sub_surface_b, s_vals, colors[0])

    # Add the control points polygon for the original surface.
    nodes = surface._nodes[(0, 2, 4, 5, 0), :]
    add_patch(ax, nodes, colors[2], with_nodes=False)

    # Add the control points polygon for the sub-surface.
    nodes = sub_surface_b._nodes[(0, 1, 2, 5, 3, 0), :]
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


def surface_locate(surface, point):
    """Image for :meth`.Surface.locate` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(256)
    ax.plot(point[:, 0], point[:, 1], color='black',
            linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.1875, 1.0625)
    save_image(ax.figure, 'surface_locate.png')


def _pretty_value(float_val):
    as_frac = fractions.Fraction(float_val)
    if as_frac.denominator == 1:
        return as_frac.numerator

    frac_tex = r'\frac{{{}}}{{{}}}'.format(
        abs(as_frac.numerator), as_frac.denominator)
    if float_val < 0.0:
        frac_tex = '-' + frac_tex
    if len(frac_tex) < 20:
        return frac_tex
    else:
        return '{:g}'.format(float_val)


def curve_specialize(curve, new_curve):
    """Image for :meth`.Curve.specialize` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256)
    interval = r'$\left[{}, {}\right]$'.format(
        _pretty_value(curve.start), _pretty_value(curve.end))
    line = ax.lines[-1]
    line.set_label(interval)
    color1 = line.get_color()

    new_curve.plot(256, ax=ax)
    interval = r'$\left[{}, {}\right]$'.format(
        _pretty_value(new_curve.start), _pretty_value(new_curve.end))
    line = ax.lines[-1]
    line.set_label(interval)

    ax.plot(curve._nodes[(0, -1), 0], curve._nodes[(0, -1), 1],
            color=color1, linestyle='None', marker='o')
    ax.plot(new_curve._nodes[(0, -1), 0], new_curve._nodes[(0, -1), 1],
            color=line.get_color(), linestyle='None', marker='o')

    ax.legend(loc='lower right', fontsize=12)
    ax.axis('scaled')
    ax.set_xlim(-0.375, 1.125)
    ax.set_ylim(-0.75, 0.625)
    save_image(ax.figure, 'curve_specialize.png')


def surface_width1(surface):
    """Image for :meth`.Surface.width` docstring."""
    if NO_IMAGES:
        return

    ax = surface.plot(2)
    ax.plot([surface.base_x], [surface.base_y],
            color='black', linestyle='None', marker='o')

    arrow_props = {
        'arrowstyle': ']-',
        'linewidth': 1,
    }
    ax.annotate(
        '', xy=(0.0, -0.09375), xycoords='data',
        xytext=(1.0, -0.09375), textcoords='data',
        color='black', arrowprops=arrow_props)

    ax.text(0.5, -0.125, '$1$', fontsize=20,
            verticalalignment='top', horizontalalignment='center')

    ax.axis('scaled')
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.25, 1.0625)
    save_image(ax.figure, 'surface_width1.png')


def surface_width2(sub_surface_b, sub_surface_c):
    """Image for :meth`.Surface.width` docstring."""
    if NO_IMAGES:
        return

    arrow_props = {
        'arrowstyle': ']-',
        'linewidth': 1,
    }

    figure, (ax1, ax2) = plt.subplots(1, 2)
    for ax in (ax1, ax2):
        ax.plot([1, 0, 0], [0, 0, 1],
                color='black', linestyle='dashed', alpha=0.5)
    # Full hypotenuse on RHS plot
    ax2.plot([0, 1], [1, 0],
             color='black', linestyle='dashed', alpha=0.5)
    # Leave space in hypotenuse for text on LHS plot
    ax1.plot([0, 0.21875], [1, 0.78125],
             color='black', linestyle='dashed', alpha=0.5)
    ax1.plot([0.21875, 0.3125], [0.78125, 0.6875],
             color='black', linestyle='dashed', alpha=0.1875)
    ax1.plot([0.3125, 1], [0.6875, 0],
             color='black', linestyle='dashed', alpha=0.5)

    sub_surface_b.plot(2, ax=ax1)
    ax1.plot([sub_surface_b.base_x], [sub_surface_b.base_y],
             color='black', linestyle='None', marker='o')
    ax1.annotate(
        '', xy=(0.5, 0.59375), xycoords='data',
        xytext=(0.0, 0.59375), textcoords='data',
        color='black', arrowprops=arrow_props)
    ax1.text(0.21875, 0.78125, r'$-\frac{1}{2}$', fontsize=20,
             verticalalignment='top', horizontalalignment='center')

    sub_surface_c.plot(2, ax=ax2)
    ax2.plot([sub_surface_c.base_x], [sub_surface_c.base_y],
             color='black', linestyle='None', marker='o')
    ax2.annotate(
        '', xy=(0.5, -0.09375), xycoords='data',
        xytext=(1.0, -0.09375), textcoords='data',
        color='black', arrowprops=arrow_props)
    ax2.text(0.75, -0.125, r'$\frac{1}{2}$', fontsize=20,
             verticalalignment='top', horizontalalignment='center')

    for ax in (ax1, ax2):
        ax.axis('scaled')
        ax.set_xlim(-0.0625, 1.0625)
        ax.set_ylim(-0.3125, 1.0625)
        save_image(ax.figure, 'surface_width2.png')


def newton_refine_surface(surface, x_val, y_val, s, t, new_s, new_t):
    """Image for :func:`._surface_helpers.newton_refine` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2)

    # Plot features of the parameter space in ax1.
    tri_surf = bezier.Surface.from_nodes(np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ]))
    tri_surf.plot(2, ax=ax1)
    ax1.plot([0.25], [0.5], marker='H')
    ax1.plot([s], [t], color='black',
             linestyle='None', marker='o')
    ax1.plot([new_s], [new_t],
             color='black', linestyle='None', marker='o',
             markeredgewidth=1, markerfacecolor='None')

    # Plot the equivalent output in ax2.
    surface.plot(256, ax=ax2)
    points = surface.evaluate_multi(np.array([
        [s, t],
        [new_s, new_t],
    ]))
    ax2.plot([x_val], [y_val], marker='H')
    ax2.plot(points[[0], 0], points[[0], 1],
             color='black', linestyle='None', marker='o')
    ax2.plot(points[[1], 0], points[[1], 1],
             color='black', linestyle='None', marker='o',
             markeredgewidth=1, markerfacecolor='None')

    # Set the axis bounds / scaling.
    ax1.axis('scaled')
    ax1.set_xlim(-0.0625, 1.0625)
    ax1.set_ylim(-0.0625, 1.0625)

    ax2.axis('scaled')
    ax2.set_xlim(-0.125, 2.125)
    ax2.set_ylim(-0.125, 2.125)

    save_image(figure, 'newton_refine_surface.png')


def classify_help(s, curve1, surface1, curve2, surface2, interior, ax=None):
    assert surface1.is_valid
    edge1, _, _ = surface1.edges
    assert np.all(edge1._nodes == curve1._nodes)

    assert surface2.is_valid
    edge2, _, _ = surface2.edges
    assert np.all(edge2._nodes == curve2._nodes)

    ax = surface1.plot(256, ax=ax)
    # Manually reduce the alpha on the surface patch(es) and
    # remove the lines we aren't using.
    ax.patches[-1].set_alpha(0.1875)
    color1 = ax.lines[-1].get_color()
    ax.lines[-1].remove()
    ax.lines[-1].remove()

    surface2.plot(256, ax=ax)
    ax.patches[-1].set_alpha(0.1875)
    color2 = ax.lines[-1].get_color()
    ax.lines[-1].remove()
    ax.lines[-1].remove()

    (int_x, int_y), = curve1.evaluate(s)
    if interior == 0:
        color = color1
    elif interior == 1:
        color = color2
    else:
        color = None
    ax.plot([int_x], [int_y],
            color=color, linestyle='None', marker='o')

    ax.axis('scaled')
    return ax


def classify_intersection1(s, curve1, tangent1, curve2, tangent2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [1.0, 0.0],
        [1.75, 0.25],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.5, 1.5],
        [1.0, 2.0],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [0.0, 0.0],
        [1.6875, 0.0625],
        [2.0, 0.5],
        [0.25, 1.0],
        [1.25, 1.25],
        [0.5, 2.0],
    ]))

    ax = classify_help(s, curve1, surface1, curve2, surface2, 0)
    (int_x, int_y), = curve1.evaluate(s)

    # Remove the alpha from the color
    color1 = ax.patches[0].get_facecolor()[:3]
    color2 = ax.patches[1].get_facecolor()[:3]
    ax.plot([int_x, int_x + tangent1[0, 0]],
            [int_y, int_y + tangent1[0, 1]],
            color=color1, linestyle='dashed')
    ax.plot([int_x, int_x + tangent2[0, 0]],
            [int_y, int_y + tangent2[0, 1]],
            color=color2, linestyle='dashed')
    ax.plot([int_x], [int_y],
            color=color1, linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 2.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, 'classify_intersection1.png')


def classify_intersection2(s, curve1, curve2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [1.0, 0.0],
        [1.5, 1.0],
        [2.0, 0.0],
        [1.25, 1.0],
        [1.75, 1.0],
        [1.5, 2.0],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [0.0, 0.0],
        [1.5, 1.0],
        [3.0, 0.0],
        [0.75, 2.0],
        [2.25, 2.0],
        [1.5, 4.0],
    ]))

    ax = classify_help(s, curve1, surface1, curve2, surface2, 1)
    ax.set_xlim(-0.0625, 3.0625)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'classify_intersection2.png')


def classify_intersection3(s, curve1, curve2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [2.0, 0.0],
        [1.5, 1.0],
        [1.0, 0.0],
        [1.75, -1.0],
        [1.25, -1.0],
        [1.5, -2.0],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [3.0, 0.0],
        [1.5, 1.0],
        [0.0, 0.0],
        [2.25, -2.0],
        [0.75, -2.0],
        [1.5, -4.0],
    ]))

    ax = classify_help(s, curve1, surface1, curve2, surface2, 0)
    ax.set_xlim(-0.0625, 3.0625)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'classify_intersection3.png')


def classify_intersection4(s, curve1, curve2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [2.0, 0.0],
        [1.5, 1.0],
        [1.0, 0.0],
        [1.75, -1.0],
        [1.25, -1.0],
        [1.5, -2.0],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [0.0, 0.0],
        [1.5, 1.0],
        [3.0, 0.0],
        [0.75, 2.0],
        [2.25, 2.0],
        [1.5, 4.0],
    ]))

    ax = classify_help(s, curve1, surface1, curve2, surface2, None)
    ax.set_xlim(-0.0625, 3.0625)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'classify_intersection4.png')


def classify_intersection5(s, curve1, curve2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [1.0, 0.0],
        [1.5, 1.0],
        [2.0, 0.0],
        [1.25, 0.9375],
        [1.75, 0.9375],
        [1.5, 1.875],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [3.0, 0.0],
        [1.5, 1.0],
        [0.0, 0.0],
        [2.25, -2.0],
        [0.75, -2.0],
        [1.5, -4.0],
    ]))

    figure, (ax1, ax2) = plt.subplots(2, 1)
    classify_help(s, curve1, surface1, curve2, surface2, 0, ax=ax1)
    classify_help(s, curve1, surface1, curve2, surface2, 1, ax=ax2)

    # Remove the alpha from the color
    color1 = ax1.patches[0].get_facecolor()[:3]
    color2 = ax1.patches[1].get_facecolor()[:3]

    # Now add the "degenerate" intersection polygons. The first
    # comes from specializing to
    # left1(0.5, 1.0)-left2(0.0, 0.25)-right1(0.375, 0.5)
    surface3 = bezier.Surface.from_nodes(np.array([
        [1.5, 0.5],
        [1.75, 0.5],
        [2.0, 0.0],
        [1.6875, 0.5],
        [1.9375, 0.234375],
        [1.875, 0.46875],
    ]))
    # NOTE: We don't require the intersection polygon be valid.
    surface3.plot(256, ax=ax1)

    # The second comes from specializing to
    # left1(0.0, 0.5)-right1(0.5, 0.625)-left3(0.75, 1.0)
    surface4 = bezier.Surface.from_nodes(np.array([
        [1.0, 0.0],
        [1.25, 0.5],
        [1.5, 0.5],
        [1.0625, 0.234375],
        [1.3125, 0.5],
        [1.125, 0.46875],
    ]))
    # NOTE: We don't require the intersection polygon be valid.
    surface4.plot(256, ax=ax2)

    (int_x, int_y), = curve1.evaluate(s)
    ax1.plot([int_x], [int_y],
             color=color1, linestyle='None', marker='o')
    ax2.plot([int_x], [int_y],
             color=color2, linestyle='None', marker='o')

    for ax in (ax1, ax2):
        ax.axis('scaled')
        ax.set_xlim(-0.0625, 3.0625)
        ax.set_ylim(-0.0625, 0.5625)

    plt.setp(ax1.get_xticklabels(), visible=False)
    figure.tight_layout(h_pad=-7.0)
    save_image(figure, 'classify_intersection5.png')


def classify_intersection6(s, curve1, curve2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [0.375, 0.0625],
        [-0.125, -0.0625],
        [-0.125, 0.0625],
        [0.1875, 0.15625],
        [-0.0625, 0.15625],
        [0.0, 0.25],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [0.75, 0.25],
        [-0.25, -0.25],
        [-0.25, 0.25],
        [0.625, 0.625],
        [0.125, 0.625],
        [0.5, 1.0],
    ]))

    ax = classify_help(s, curve1, surface1, curve2, surface2, None)
    ax.set_xlim(-0.3125, 1.0625)
    ax.set_ylim(-0.0625, 0.3125)
    save_image(ax.figure, 'classify_intersection6.png')


def classify_intersection7(s, curve1a, curve1b, curve2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    surface1 = bezier.Surface.from_nodes(np.array([
        [0.0, 0.0],
        [4.5, 0.0],
        [9.0, 2.25],
        [0.0, 1.25],
        [4.5, 2.375],
        [0.0, 2.5],
    ]))
    surface2 = bezier.Surface.from_nodes(np.array([
        [11.25, 0.0],
        [9.0, 4.5],
        [2.75, 1.0],
        [8.125, -0.75],
        [3.875, -0.25],
        [5.0, -1.5],
    ]))

    figure, (ax1, ax2) = plt.subplots(2, 1)
    classify_help(s, curve1a, surface1, curve2, surface2, None, ax=ax1)
    surface1._nodes = surface1._nodes[(2, 4, 5, 1, 3, 0), :]
    surface1._edges = None
    classify_help(0.0, curve1b, surface1, curve2, surface2, 0, ax=ax2)

    for ax in (ax1, ax2):
        ax.set_xlim(-0.125, 11.5)
        ax.set_ylim(-0.125, 2.625)

    plt.setp(ax1.get_xticklabels(), visible=False)
    figure.tight_layout(h_pad=-5.0)
    save_image(figure, 'classify_intersection7.png')


def get_curvature(nodes, s, tangent_vec, curvature):
    """Image for :func:`get_curvature` docstring."""
    if NO_IMAGES:
        return

    curve = bezier.Curve.from_nodes(nodes)

    # Find the center of the circle along the direction
    # perpendicular to the tangent vector (90 degree left turn).
    radius_dir = np.array([[-tangent_vec[0, 1], tangent_vec[0, 0]]])
    radius_dir /= np.linalg.norm(radius_dir, ord=2)
    point = curve.evaluate(s)
    circle_center = point + radius_dir / curvature

    # Add the curve.
    ax = curve.plot(256)
    # Add the circle.
    circle_center = circle_center.flatten()
    circle = plt.Circle(circle_center, 1.0 / abs(curvature), alpha=0.25)
    ax.add_artist(circle)
    # Add the point.
    ax.plot(point[:, 0], point[:, 1],
            color='black', marker='o', linestyle='None')

    ax.axis('scaled')
    ax.set_xlim(-0.0625, 1.0625)
    ax.set_ylim(-0.0625, 0.625)
    save_image(ax.figure, 'get_curvature.png')


def curve_locate(curve, point1, point2):
    """Image for :meth`.Curve.locate` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256)
    points = np.vstack([point1, point2])
    ax.plot(points[:, 0], points[:, 1], color='black',
            linestyle='None', marker='o')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 4.125)
    ax.set_ylim(-0.125, 1.25)
    save_image(ax.figure, 'curve_locate.png')


def newton_refine_curve(curve, point, s, new_s):
    """Image for :func:`._curve_helpers.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256)
    ax.plot(point[:, 0], point[:, 1], marker='H')
    wrong_points = curve.evaluate_multi(np.array([s, new_s]))
    ax.plot(wrong_points[[0], 0], wrong_points[[0], 1],
            color='black', linestyle='None', marker='o')
    ax.plot(wrong_points[[1], 0], wrong_points[[1], 1],
            color='black', linestyle='None', marker='o',
            markeredgewidth=1, markerfacecolor='None')

    # Set the axis bounds / scaling.
    ax.axis('scaled')
    ax.set_xlim(-0.125, 3.125)
    ax.set_ylim(-0.125, 1.375)

    save_image(ax.figure, 'newton_refine_curve.png')


def newton_refine_curve_cusp(curve, s_vals):
    """Image for :func:`._curve_helpers.newton_refine` docstring."""
    if NO_IMAGES:
        return

    ax = curve.plot(256)
    ax.lines[-1].zorder = 1

    points = curve.evaluate_multi(np.array(s_vals))
    colors = seaborn.dark_palette('blue', 6)
    ax.scatter(points[:, 0], points[:, 1], c=colors,
               s=20, alpha=0.75, zorder=2)

    # Set the axis bounds / scaling.
    ax.axis('scaled')
    ax.set_xlim(-0.125, 6.125)
    ax.set_ylim(-3.125, 3.125)

    save_image(ax.figure, 'newton_refine_curve_cusp.png')


def classify_intersection8(s, curve1, surface1, curve2, surface2):
    """Image for :func:`._surface_helpers.classify_intersection` docstring."""
    if NO_IMAGES:
        return

    ax = classify_help(s, curve1, surface1, curve2, surface2, None)
    ax.set_xlim(-1.125, 1.125)
    ax.set_ylim(-0.125, 1.125)
    save_image(ax.figure, 'classify_intersection8.png')


def curve_elevate(curve, elevated):
    """Image for :meth:`.curve.Curve.elevate` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2)

    curve.plot(256, ax=ax1)
    color = ax1.lines[-1].get_color()
    add_patch(ax1, curve._nodes, color)

    elevated.plot(256, ax=ax2)
    color = ax2.lines[-1].get_color()
    add_patch(ax2, elevated._nodes, color)

    ax1.axis('scaled')
    ax2.axis('scaled')
    _plot_helpers.add_plot_boundary(ax1)
    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*ax1.get_ylim())

    save_image(figure, 'curve_elevate.png')


def surface_elevate(surface, elevated):
    """Image for :meth:`.surface.Surface.elevate` docstring."""
    if NO_IMAGES:
        return

    figure, (ax1, ax2) = plt.subplots(1, 2)

    surface.plot(256, ax=ax1)
    color = ax1.lines[-1].get_color()
    nodes = surface._nodes[(0, 1, 2, 4, 5), :]
    add_patch(ax1, nodes, color)

    elevated.plot(256, ax=ax2)
    color = ax2.lines[-1].get_color()
    nodes = elevated._nodes[(0, 1, 2, 3, 6, 8, 9), :]
    add_patch(ax2, nodes, color)

    ax1.axis('scaled')
    ax2.axis('scaled')
    _plot_helpers.add_plot_boundary(ax1)
    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*ax1.get_ylim())

    save_image(figure, 'surface_elevate.png')


def unit_triangle():
    """Image for :class:`.surface.Surface` docstring."""
    if NO_IMAGES:
        return

    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ])
    surface = bezier.Surface(nodes, degree=1)

    ax = surface.plot(256)

    ax.axis('scaled')
    _plot_helpers.add_plot_boundary(ax)

    save_image(ax.figure, 'unit_triangle.png')
