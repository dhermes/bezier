# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

import itertools

import matplotlib.pyplot as plt
import numpy as np
import pytest
try:
    import seaborn  # pylint: disable=unused-import
except ImportError:
    pass
import six

import bezier
from bezier import _intersection_helpers

import runtime_utils


CONFIG = runtime_utils.Config()

# F1L = sympy.Matrix([[s, t]])
SURFACE1L = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
]))
# F2L = sympy.Matrix([[
#     (9 * s + 2 * t - 1) / 8,
#     (9 * s + 7 * t - 1) / 16,
# ]])
SURFACE2L = bezier.Surface(np.array([
    [-0.125, -0.0625],
    [1.0, 0.5],
    [0.125, 0.375],
]))
# F3L = sympy.Matrix([[(4 * s + 1) / 4, (8 * t - 1) / 8]])
SURFACE3L = bezier.Surface(np.array([
    [0.25, -0.125],
    [1.25, -0.125],
    [0.25, 0.875],
]))
# F4L = sympy.Matrix([[5 * (2 * s + t - 1) / 4, 35 * t / 16]])
SURFACE4L = bezier.Surface(np.array([
    [-1.25, 0.0],
    [1.25, 0.0],
    [0.0, 2.1875],
]))
# F5L = sympy.Matrix([[
#     (21 * s - 21 * t + 4) / 8,
#     (21 * s + 21 * t - 5) / 32,
# ]])
SURFACE5L = bezier.Surface(np.array([
    [0.5, -0.15625],
    [3.125, 0.5],
    [-2.125, 0.5],
]))
# F6L = sympy.Matrix([[(1 - 4 * s) / 4, (1 - 4 * t) / 4]])
SURFACE6L = bezier.Surface(np.array([
    [0.25, 0.25],
    [-0.75, 0.25],
    [0.25, -0.75],
]))

# F1Q = sympy.Matrix([[
#     (2 * s - t**2 + t) / 2,
#     (s**2 + 2 * s * t - s + 2 * t) / 2,
# ]])
SURFACE1Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, -0.25],
    [1.0, 0.0],
    [0.25, 0.5],
    [0.75, 0.75],
    [0.0, 1.0],
]))
# F2Q = sympy.Matrix([[
#     (3 * s**2 + 6 * s * t + 5 * s + 8 * t - 2) / 8,
#     (3 * s**2 - 3 * t**2 - 11 * s + 3 * t + 6) / 8,
# ]])
SURFACE2Q = bezier.Surface(np.array([
    [-0.25, 0.75],
    [0.0625, 0.0625],
    [0.75, -0.25],
    [0.25, 0.9375],
    [0.9375, 0.25],
    [0.75, 0.75],
]))
# F3Q = sympy.Matrix([[
#     (11 * s * t + 5 * t**2 + 8 * s + 3 * t) / 8,
#     -(4 * s**2 + 7 * s * t + 7 * t**2 - 4 * s - 15 * t) / 8,
# ]])
SURFACE3Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.25],
    [1.0, 0.0],
    [0.1875, 0.9375],
    [1.375, 0.75],
    [1.0, 1.0],
]))
# F4Q = sympy.Matrix([[
#     -(2 * s * t + t**2 - 16 * s - 4 * t) / 16,
#     (2 * s**2 + 2 * s * t - 2 * s + 3 * t + 1) / 4,
# ]])
SURFACE4Q = bezier.Surface(np.array([
    [0.0, 0.25],
    [0.5, 0.0],
    [1.0, 0.25],
    [0.125, 0.625],
    [0.5625, 0.625],
    [0.1875, 1.0],
]))
# F5Q = sympy.Matrix([[
#     -(s**2 + s * t - 8 * s - 3 * t - 2) / 8,
#     -(25 * s**2 + 20 * s * t - t**2 - 34 * s - 28 * t - 3) / 32,
# ]])
SURFACE5Q = bezier.Surface(np.array([
    [0.25, 0.09375],
    [0.75, 0.625],
    [1.125, 0.375],
    [0.4375, 0.53125],
    [0.875, 0.75],
    [0.625, 1.0],
]))
# F6Q = sympy.Matrix([[
#     (s**2 + s * t + 3 * s + 2 * t) / 4,
#     -(3 * s**2 + 3 * s * t - 3 * s - 4 * t) / 4,
# ]])
SURFACE6Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.375, 0.375],
    [1.0, 0.0],
    [0.25, 0.5],
    [0.75, 0.5],
    [0.5, 1.0],
]))
# F7Q = sympy.Matrix([[
#     -(s**2 + s * t + 3 * s + 2 * t - 4) / 4,
#     (8 * s**2 + 8 * s * t - 8 * s - 9 * t + 3) / 8,
# ]])
SURFACE7Q = bezier.Surface(np.array([
    [1.0, 0.375],
    [0.625, -0.125],
    [0.0, 0.375],
    [0.75, -0.1875],
    [0.25, -0.1875],
    [0.5, -0.75],
]))
# F8Q = sympy.Matrix([[
#     s * (7 * t + 1),
#     t * (7 * s + 1),
# ]])
SURFACE8Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.0],
    [1.0, 0.0],
    [0.0, 0.5],
    [4.0, 4.0],
    [0.0, 1.0],
]))
# F9Q = sympy.Matrix([[
#     (14 * s * t + 14 * t**2 + 2 * s - 12 * t + 3) / 2,
#     -t * (7 * s + 7 * t - 8),
# ]])
SURFACE9Q = bezier.Surface(np.array([
    [1.5, 0.0],
    [2.0, 0.0],
    [2.5, 0.0],
    [-1.5, 4.0],
    [2.5, 0.5],
    [2.5, 1.0],
]))
# F10Q = sympy.Matrix([[
#     -(2 * s + t - 2) / 2,
#     -(2 * s**2 + 2 * s*t - 2 * s + 3 * t) / 4,
# ]])
SURFACE10Q = bezier.Surface(np.array([
    [1.0, 0.0],
    [0.5, 0.25],
    [0.0, 0.0],
    [0.75, -0.375],
    [0.25, -0.375],
    [0.5, -0.75],
]))
# F11Q = sympy.Matrix([[
#     (16 * s - t**2 + 4 * t) / 16,
#     (8 * s**2 + 8 * s * t - 8 * s + 12 * t + 3) / 16,
# ]])
SURFACE11Q = bezier.Surface(np.array([
    [0.0, 0.1875],
    [0.5, -0.0625],
    [1.0, 0.1875],
    [0.125, 0.5625],
    [0.625, 0.5625],
    [0.1875, 0.9375],
]))
# F12Q = sympy.Matrix([[
#     -(2 * s + t - 2) / 2,
#     -(8 * s**2 + 8 * s * t - 8 * s + 12 * t - 1) / 16,
# ]])
SURFACE12Q = bezier.Surface(np.array([
    [1.0, 0.0625],
    [0.5, 0.3125],
    [0.0, 0.0625],
    [0.75, -0.3125],
    [0.25, -0.3125],
    [0.5, -0.6875],
]))
# F13Q = sympy.Matrix([[
#     (2 * s + t + 1) / 4,
#     (4 * s**2 + 4 * s * t + t**2 - 4 * s + 14 * t + 5) / 32,
# ]])
# NOTE: The bottom edge of this surface lies entirely on the
#       bottom edge of SURFACE4Q.
SURFACE13Q = bezier.Surface(np.array([
    [0.25, 0.15625],
    [0.5, 0.09375],
    [0.75, 0.15625],
    [0.375, 0.375],
    [0.625, 0.375],
    [0.5, 0.625],
]))
# F14Q = F13Q + sympy.Matrix([[0, 1]]) / 16
SURFACE14Q = bezier.Surface(np.array([
    [0.25, 0.21875],
    [0.5, 0.15625],
    [0.75, 0.21875],
    [0.375, 0.4375],
    [0.625, 0.4375],
    [0.5, 0.6875],
]))
# F15Q = sympy.Matrix([[
#     -2 * (s + t) * (t - 1),
#     (7 * s + 8 * t) / 4,
# ]])
SURFACE15Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 0.875],
    [2.0, 1.75],
    [1.0, 1.0],
    [1.0, 1.875],
    [0.0, 2.0],
]))
# F16Q = sympy.Matrix([[
#     2 * s**2 + 2 * s * t - 2 * s - 2 * t + 1,
#     (8 * s + 7 * t) / 4,
# ]])
SURFACE16Q = bezier.Surface(np.array([
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 2.0],
    [0.0, 0.875],
    [0.0, 1.875],
    [-1.0, 1.75],
]))
# F17Q = sympy.Matrix([[
#     (19 * s - 19 * t + 32) / 64,
#     (2 *s * t + 5 * s + 5 * t - 6) / 8,
# ]])
SURFACE17Q = bezier.Surface(np.array([
    [0.5, -0.75],
    [0.6484375, -0.4375],
    [0.796875, -0.125],
    [0.3515625, -0.4375],
    [0.5, 0.0],
    [0.203125, -0.125],
]))
# F18Q = sympy.Matrix([[
#     -(4 * s + t - 3) / 4,
#     -(16 * s**2 + 16 * s * t - 8 * s + 27 * t - 3) / 32,
# ]])
SURFACE18Q = bezier.Surface(np.array([
    [0.75, 0.09375],
    [0.25, 0.21875],
    [-0.25, -0.15625],
    [0.625, -0.328125],
    [0.125, -0.453125],
    [0.5, -0.75],
]))
# F19Q = sympy.Matrix([[
#     (t - s) / 2,
#     -(s**2 + s * t + 2 * s + 6 * t) / 8,
# ]])
SURFACE19Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [-0.25, -0.125],
    [-0.5, -0.375],
    [0.25, -0.375],
    [0.0, -0.5625],
    [0.5, -0.75],
]))
# F20Q = sympy.Matrix([[
#     -s**2 - s * t - 2 * t + 1, -s * (s + t - 2),
# ]])
SURFACE20Q = bezier.Surface(np.array([
    [1.0, 0.0],
    [1.0, 1.0],
    [0.0, 1.0],
    [0.0, 0.0],
    [-0.5, 0.5],
    [-1.0, 0.0],
]))
# F21Q = sympy.Matrix([[
#     -(3 * s**2 + 3 * s * t + 3 * t - 2) / 2,
#     -(6 * s**2 + 6 * s * t - 12 * s - t) / 4,
# ]])
SURFACE21Q = bezier.Surface(np.array([
    [1.0, 0.0],
    [1.0, 1.5],
    [-0.5, 1.5],
    [0.25, 0.125],
    [-0.5, 0.875],
    [-0.5, 0.25],
]))
# F22Q = sympy.Matrix([[
#     -5 * (2 * t - 5) * (2 * s + t - 1) / 16,
#     -5 * (16 * s**2 + 16 * s * t + 2 * t**2 - 16 * s - 37 * t + 4) / 64,
# ]])
SURFACE22Q = bezier.Surface(np.array([
    [-1.5625, -0.3125],
    [0.0, 0.3125],
    [1.5625, -0.3125],
    [-0.46875, 1.1328125],
    [0.46875, 1.1328125],
    [0.0, 2.421875],
]))
# F23Q = sympy.Matrix([[
#     (t + 2) * (2 * s + t - 1) / 2,
#     (4 * s**2 + 4 * s * t - 3 * t**2 - 4 * s + 17 * t + 1) / 8,
# ]])
SURFACE23Q = bezier.Surface(np.array([
    [-1.0, 0.125],
    [0.0, -0.125],
    [1.0, 0.125],
    [-0.75, 1.1875],
    [0.75, 1.1875],
    [0.0, 1.875],
]))
# F24Q = sympy.Matrix([[
#     -(42 * s**2 - 10 * s * t - 33 * t**2 + 16 * s + 128 * t - 96) / 128,
#     (3 * s**2 + 24 * s * t + 11 * t**2 - 18 * t + 25) / 32,
# ]])
SURFACE24Q = bezier.Surface(np.array([
    [0.75, 0.78125],
    [0.6875, 0.78125],
    [0.296875, 0.875],
    [0.25, 0.5],
    [0.2265625, 0.875],
    [0.0078125, 0.5625],
]))
# F25Q = sympy.Matrix([[
#     -(s**2 + 38 * s * t + 44 * t**2 + 40 * s + 8 * t - 62) / 64,
#     (2 * s**2 + 12 * s * t + t**2 - 24 * s - 56 * t + 62) / 64,
# ]])
SURFACE25Q = bezier.Surface(np.array([
    [0.96875, 0.96875],
    [0.65625, 0.78125],
    [0.328125, 0.625],
    [0.90625, 0.53125],
    [0.296875, 0.4375],
    [0.15625, 0.109375],
]))
# F26Q = sympy.Matrix([[s * (t + 1), t * (s + 1)]])
SURFACE26Q = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.0],
    [1.0, 0.0],
    [0.0, 0.5],
    [1.0, 1.0],
    [0.0, 1.0],
]))
# F27Q = sympy.Matrix([[s * t + t**2 - 2 * t + 2, s * t + t**2 + s]])
SURFACE27Q = bezier.Surface(np.array([
    [2.0, 0.0],
    [2.0, 0.5],
    [2.0, 1.0],
    [1.0, 0.0],
    [1.5, 1.0],
    [1.0, 1.0],
]))
# F28Q = sympy.Matrix([[
#     -(4 * s**2 + 6 * s * t - 3 * t**2 - 2 * s + 24 * t - 24) / 16,
#     (3 * s * t + 3 * t**2 + 8 * s - 3 * t) / 8,
# ]])
SURFACE28Q = bezier.Surface(np.array([
    [1.5, 0.0],
    [1.5625, 0.5],
    [1.375, 1.0],
    [0.75, -0.1875],
    [0.625, 0.5],
    [0.1875, 0.0],
]))
# F29Q = sympy.Matrix([[
#     (s**2 - 7 * s * t + 13 * s + 4 * t - 4) / 8,
#     -(7 * s * t - t**2 - 4 * s - 13 * t + 4) / 8,
# ]])
SURFACE29Q = bezier.Surface(np.array([
    [-0.5, -0.5],
    [0.3125, -0.25],
    [1.25, 0.0],
    [-0.25, 0.3125],
    [0.125, 0.125],
    [0.0, 1.25],
]))


def edge_pairs(surface1, surface2):
    edges1 = six.moves.map(
        _intersection_helpers.Linearization.from_shape,
        surface1.edges)
    edges2 = six.moves.map(
        _intersection_helpers.Linearization.from_shape,
        surface2.edges)
    return itertools.product(edges1, edges2)


def make_plots(surface1, surface2):
    if not CONFIG.running:
        return

    ax = surface1.plot(64)
    surface2.plot(64, ax=ax)
    plt.axis('scaled')

    if CONFIG.save_plot:
        CONFIG.save_fig()
    else:
        plt.title(CONFIG.current_test)
        plt.show()

    plt.close(ax.figure)


def surface_surface_check(surface1, surface2):
    assert surface1.is_valid
    assert surface2.is_valid

    _intersection_helpers.all_intersections(
        edge_pairs(surface1, surface2))

    make_plots(surface1, surface2)


def test_surfaces1Q_and_3Q():
    surface_surface_check(SURFACE1Q, SURFACE3Q)


def test_surfaces1L_and_3L():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE1L, SURFACE3L)

    make_plots(SURFACE1L, SURFACE3L)


def test_surfaces1Q_and_2Q():
    surface_surface_check(SURFACE1Q, SURFACE2Q)


def test_surfaces10Q_and_18Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE10Q, SURFACE18Q)

    make_plots(SURFACE10Q, SURFACE18Q)


def test_surfaces10Q_and_19Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE10Q, SURFACE19Q)

    make_plots(SURFACE10Q, SURFACE19Q)


def test_surfaces3Q_and_4Q():
    surface_surface_check(SURFACE3Q, SURFACE4Q)


def test_surfaces1Q_and_5L():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE1Q, SURFACE5L)

    make_plots(SURFACE1Q, SURFACE5L)


def test_surfaces3Q_and_5Q():
    surface_surface_check(SURFACE3Q, SURFACE5Q)


def test_surfaces1L_and_2L():
    surface_surface_check(SURFACE1L, SURFACE2L)


def test_surfaces20Q_and_21Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE20Q, SURFACE21Q)

    make_plots(SURFACE20Q, SURFACE21Q)


def test_surfaces4L_and_22Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE4L, SURFACE22Q)

    make_plots(SURFACE4L, SURFACE22Q)


def test_surfaces4L_and_23Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE4L, SURFACE23Q)

    make_plots(SURFACE4L, SURFACE23Q)


def test_surfaces6Q_and_7Q():
    surface_surface_check(SURFACE6Q, SURFACE7Q)


def test_surfaces8Q_and_9Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE8Q, SURFACE9Q)

    make_plots(SURFACE8Q, SURFACE9Q)


def test_surfaces4Q_and_10Q():
    surface_surface_check(SURFACE4Q, SURFACE10Q)


def test_surfaces11Q_and_12Q():
    surface_surface_check(SURFACE11Q, SURFACE12Q)


def test_surfaces3Q_and_13Q():
    surface_surface_check(SURFACE3Q, SURFACE13Q)


def test_surfaces10Q_and_17Q():
    surface_surface_check(SURFACE10Q, SURFACE17Q)


def test_surfaces3Q_and_14Q():
    surface_surface_check(SURFACE3Q, SURFACE14Q)


def test_surfaces15Q_and_16Q():
    surface_surface_check(SURFACE15Q, SURFACE16Q)


def test_surfaces24Q_and_25Q():
    surface_surface_check(SURFACE24Q, SURFACE25Q)


def test_surfaces1L_and_6L():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE1L, SURFACE6L)

    make_plots(SURFACE1L, SURFACE6L)


def test_surfaces26Q_and_27Q():
    with pytest.raises(NotImplementedError):
        surface_surface_check(SURFACE26Q, SURFACE27Q)

    make_plots(SURFACE26Q, SURFACE27Q)


def test_surfaces1L_and_28Q():
    surface_surface_check(SURFACE1L, SURFACE28Q)


def test_surfaces1L_and_29Q():
    surface_surface_check(SURFACE1L, SURFACE29Q)


if __name__ == '__main__':
    CONFIG.run(globals())
