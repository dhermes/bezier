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


def save_image(figure, filename):
    """Save an image to the docs images directory.

    Args:
        filename (str): The name of the file (not containing
            directory info).
    """
    path = os.path.join(IMAGES_DIR, filename)
    figure.savefig(path, bbox_inches='tight')
    print('Saved {}'.format(filename))
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
               s=30, alpha=0.75, zorder=2)

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


def add_patch(ax, nodes, color):
    path = _path_mod.Path(nodes)
    patch = patches.PathPatch(
        path, facecolor=color, alpha=0.6)
    ax.add_patch(patch)
    ax.plot(nodes[:, 0], nodes[:, 1], color='black',
            linestyle='None', marker='o')


def curve_constructor():
    """Image for :class`.Curve` docstring."""
    nodes = np.array([
        [0.0  , 0.0],
        [0.625, 0.5],
        [1.0  , 0.5],
    ])
    curve = bezier.Curve(nodes)

    ax = curve.plot(256)
    line = ax.lines[0]
    ax.plot(nodes[:, 0], nodes[:, 1], linestyle='None',
            marker='o', color='black')
    add_patch(ax, nodes, line.get_color())

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'curve_constructor.png')


def curve_evaluate():
    """Image for :meth`.Curve.evaluate` docstring."""
    curve = bezier.Curve(np.array([
        [0.0  , 0.0],
        [0.625, 0.5],
        [1.0  , 0.5],
    ]))

    ax = curve.plot(256)
    points = curve.evaluate_multi(np.array([0.75]))
    ax.plot(points[:, 0], points[:, 1], linestyle='None',
            marker='o', color='black')

    ax.axis('scaled')
    ax.set_xlim(-0.125, 1.125)
    ax.set_ylim(-0.0625, 0.5625)
    save_image(ax.figure, 'curve_evaluate.png')


def curve_subdivide():
    """Image for :meth`.Curve.subdivide` docstring."""
    curve = bezier.Curve(np.array([
        [0.0, 0.0],
        [1.25, 3.0],
        [2.0, 1.0],
    ]))
    left, right = curve.subdivide()

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


def curve_intersect():
    """Image for :meth`.Curve.intersect` docstring."""
    curve1 = bezier.Curve(np.array([
        [0.0, 0.0],
        [0.375, 0.75],
        [0.75, 0.375],
    ]))
    curve2 = bezier.Curve(np.array([
        [0.5, 0.0],
        [0.5, 0.75],
    ]))

    ax = curve1.plot(256)
    curve2.plot(256, ax=ax)
    ax.plot([0.5], [0.5], linestyle='None',
            marker='o', color='black')

    ax.axis('scaled')
    ax.set_xlim(0.0, 0.75)
    ax.set_ylim(0.0, 0.75)
    save_image(ax.figure, 'curve_intersect.png')


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
    curve_constructor()
    curve_evaluate()
    curve_subdivide()
    curve_intersect()


if __name__ == '__main__':
    main()
