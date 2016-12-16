# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

import matplotlib.pyplot as plt
import numpy as np

import bezier
from bezier import _plot_helpers

import runtime_utils


CONFIG = runtime_utils.Config()


# F1 = sympy.Matrix([[s, t]])
SURFACE1 = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
]))
# F2 = sympy.Matrix([[
#     (-t**2 + 2 * s + t) / 2, (s**2 + 2 * s * t - s + 2 * t) / 2]])
SURFACE2 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, -0.25],
    [1.0, 0.0],
    [0.25, 0.5],
    [0.75, 0.75],
    [0.0, 1.0],
]))
# F3 = sympy.Matrix([[
#     -(2 * s * t - 4 * s - t) / 4, (s**2 - s * t + 4 * t) / 4]])
SURFACE3 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.0],
    [1.0, 0.25],
    [0.125, 0.5],
    [0.375, 0.375],
    [0.25, 1.0],
]))
# F4 = sympy.Matrix([[2 * (s + 2 * t) * (1 - t), 2 * t * (s + 1)]])
SURFACE4 = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [2.0, 0.0],
    [2.0, 1.0],
    [2.0, 2.0],
    [0.0, 2.0],
]))

POINTS = np.array([
    [0.0, 0.0],
    [0.25, 0.25],
    [0.59375, 0.25],
    [0.265625, 0.73046875],
    [1.25, 1.25],
])


def make_plot(surface, point):
    if not CONFIG.running:
        return

    ax = surface.plot(64)
    ax.plot(point[:, 0], point[:, 1], color='black',
            marker='o', linestyle='None')

    ax.axis('scaled')
    _plot_helpers.add_plot_boundary(ax)

    if CONFIG.save_plot:
        CONFIG.save_fig()
    else:
        plt.title(CONFIG.current_test)
        plt.show()

    plt.close(ax.figure)


def check_point(surface, point_ind, expected_s, expected_t):
    point = POINTS[[point_ind], :]
    if expected_s is None and expected_t is None:
        assert surface.locate(point) is None
    else:
        s, t = surface.locate(point)

        CONFIG.assert_close(s, expected_s)
        CONFIG.assert_close(t, expected_t)

    make_plot(surface, point)


def test_surface1_and_point0():
    check_point(SURFACE1, 0, 0.0, 0.0)


def test_surface1_and_point1():
    check_point(SURFACE1, 1, 0.25, 0.25)


def test_surface2_and_point1():
    s, _ = runtime_utils.real_roots([4, -32, -56, -24, 5])
    _, t = runtime_utils.real_roots([4, 8, -16, 44, -11])
    check_point(SURFACE2, 1, s, t)


def test_surface2_and_point2():
    check_point(SURFACE2, 2, 0.5, 0.25)


def test_surface2_and_point4():
    check_point(SURFACE2, 4, None, None)


def test_surface3_and_point1():
    s, = runtime_utils.real_roots([2, -5, 15, -3])
    t, = runtime_utils.real_roots([14, -61, 74, -15])
    check_point(SURFACE3, 1, s, t)


def test_surface3_and_point3():
    check_point(SURFACE3, 3, 0.125, 0.75)


def test_surface4_and_point2():
    s, = runtime_utils.real_roots([64, 101, 34, -5])
    t, = runtime_utils.real_roots([128, -192, 91, -8])
    check_point(SURFACE4, 2, s, t)


def test_surface4_and_point4():
    check_point(SURFACE4, 4, 0.25, 0.5)


if __name__ == '__main__':
    CONFIG.run(globals())
