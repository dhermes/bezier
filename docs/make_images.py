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


def main():
    linearization_error()
    newton_refine1()
    newton_refine2()
    newton_refine3()


if __name__ == '__main__':
    main()
