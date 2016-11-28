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
try:
    import seaborn  # pylint: disable=unused-import
except ImportError:
    pass

import bezier

import runtime_utils


CONFIG = runtime_utils.Config()
# F1 = sympy.Matrix([[
#     (2 * s - t**2 + t) / 2,
#     (s**2 + 2 * s * t - s + 2 * t) / 2,
# ]])
SURFACE1 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, -0.25],
    [1.0, 0.0],
    [0.25, 0.5],
    [0.75, 0.75],
    [0.0, 1.0],
]))
# F2 = sympy.Matrix([[
#     (3 * s**2 + 6 * s * t + 5 * s + 8 * t - 2) / 8,
#     (3 * s**2 - 3 * t**2 - 11 * s + 3 * t + 6) / 8,
# ]])
SURFACE2 = bezier.Surface(np.array([
    [-0.25, 0.75],
    [0.0625, 0.0625],
    [0.75, -0.25],
    [0.25, 0.9375],
    [0.9375, 0.25],
    [0.75, 0.75],
]))
# F3 = sympy.Matrix([[
#     (11 * s * t + 5 * t**2 + 8 * s + 3 * t) / 8,
#     -(4 * s**2 + 7 * s * t + 7 * t**2 - 4 * s - 15 * t) / 8,
# ]])
SURFACE3 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.25],
    [1.0, 0.0],
    [0.1875, 0.9375],
    [1.375, 0.75],
    [1.0, 1.0],
]))
# F4 = sympy.Matrix([[
#     -(2 * s * t + t**2 - 16 * s - 4 * t) / 16,
#     (2 * s**2 + 2 * s * t - 2 * s + 3 * t + 1) / 4,
# ]])
SURFACE4 = bezier.Surface(np.array([
    [0.0, 0.25],
    [0.5, 0.0],
    [1.0, 0.25],
    [0.125, 0.625],
    [0.5625, 0.625],
    [0.1875, 1.0],
]))
# F5 = sympy.Matrix([[
#     -(s**2 + s * t - 8 * s - 3 * t - 2) / 8,
#     -(25 * s**2 + 20 * s * t - t**2 - 34 * s - 28 * t - 3) / 32,
# ]])
SURFACE5 = bezier.Surface(np.array([
    [0.25, 0.09375],
    [0.75, 0.625],
    [1.125, 0.375],
    [0.4375, 0.53125],
    [0.875, 0.75],
    [0.625, 1.0],
]))
# F6 = sympy.Matrix([[
#     (s**2 + s * t + 3 * s + 2 * t) / 4,
#     -(3 * s**2 + 3 * s * t - 3 * s - 4 * t) / 4,
# ]])
SURFACE6 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.375, 0.375],
    [1.0, 0.0],
    [0.25, 0.5],
    [0.75, 0.5],
    [0.5, 1.0],
]))
# F7 = sympy.Matrix([[
#     -(s**2 + s * t + 3 * s + 2 * t - 4) / 4,
#     (8 * s**2 + 8 * s * t - 8 * s - 9 * t + 3) / 8,
# ]])
SURFACE7 = bezier.Surface(np.array([
    [1.0, 0.375],
    [0.625, -0.125],
    [0.0, 0.375],
    [0.75, -0.1875],
    [0.25, -0.1875],
    [0.5, -0.75],
]))
# F8 = sympy.Matrix([[
#     s * (7 * t + 1),
#     t * (7 * s + 1),
# ]])
SURFACE8 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.0],
    [1.0, 0.0],
    [0.0, 0.5],
    [4.0, 4.0],
    [0.0, 1.0],
]))
# F9 = sympy.Matrix([[
#     (14 * s * t + 14 * t**2 + 2 * s - 12 * t + 3) / 2,
#     -t * (7 * s + 7 * t - 8),
# ]])
SURFACE9 = bezier.Surface(np.array([
    [1.5, 0.0],
    [2.0, 0.0],
    [2.5, 0.0],
    [-1.5, 4.0],
    [2.5, 0.5],
    [2.5, 1.0],
]))
# F10 = sympy.Matrix([[
#     -(2 * s + t - 2) / 2,
#     -(2 * s**2 + 2 * s*t - 2 * s + 3 * t) / 4,
# ]])
SURFACE10 = bezier.Surface(np.array([
    [1.0, 0.0],
    [0.5, 0.25],
    [0.0, 0.0],
    [0.75, -0.375],
    [0.25, -0.375],
    [0.5, -0.75],
]))
# F11 = sympy.Matrix([[
#     (16 * s - t**2 + 4 * t) / 16,
#     (8 * s**2 + 8 * s * t - 8 * s + 12 * t + 3) / 16,
# ]])
SURFACE11 = bezier.Surface(np.array([
    [0.0, 0.1875],
    [0.5, -0.0625],
    [1.0, 0.1875],
    [0.125, 0.5625],
    [0.625, 0.5625],
    [0.1875, 0.9375],
]))
# F12 = sympy.Matrix([[
#     -(2 * s + t - 2) / 2,
#     -(8 * s**2 + 8 * s * t - 8 * s + 12 * t - 1) / 16,
# ]])
SURFACE12 = bezier.Surface(np.array([
    [1.0, 0.0625],
    [0.5, 0.3125],
    [0.0, 0.0625],
    [0.75, -0.3125],
    [0.25, -0.3125],
    [0.5, -0.6875],
]))
# F13 = sympy.Matrix([[
#     (2 * s + t + 1) / 4,
#     (4 * s**2 + 4 * s * t + t**2 - 4 * s + 14 * t + 5) / 32,
# ]])
# NOTE: Related to SURFACE4
SURFACE13 = bezier.Surface(np.array([
    [0.25, 0.15625],
    [0.5, 0.09375],
    [0.75, 0.15625],
    [0.375, 0.375],
    [0.625, 0.375],
    [0.5, 0.625],
]))
# F14 = F13 + sympy.Matrix([[0, 1]]) / 16
SURFACE14 = bezier.Surface(np.array([
    [0.25, 0.21875],
    [0.5, 0.15625],
    [0.75, 0.21875],
    [0.375, 0.4375],
    [0.625, 0.4375],
    [0.5, 0.6875],
]))
# F15 = sympy.Matrix([[s, t]])
SURFACE15 = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
]))
# F16 = sympy.Matrix([[
#     (9 * s + 2 * t - 1) / 8,
#     (9 * s + 7 * t - 1) / 16,
# ]])
SURFACE16 = bezier.Surface(np.array([
    [-0.125, -0.0625],
    [1.0, 0.5],
    [0.125, 0.375],
]))
# F17 = sympy.Matrix([[
#     -2 * (s + t) * (t - 1),
#     (7 * s + 8 * t) / 4,
# ]])
SURFACE17 = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 0.875],
    [2.0, 1.75],
    [1.0, 1.0],
    [1.0, 1.875],
    [0.0, 2.0],
]))
# F18 = sympy.Matrix([[
#     2 * s**2 + 2 * s * t - 2 * s - 2 * t + 1,
#     (8 * s + 7 * t) / 4,
# ]])
SURFACE18 = bezier.Surface(np.array([
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 2.0],
    [0.0, 0.875],
    [0.0, 1.875],
    [-1.0, 1.75],
]))
# F19 = sympy.Matrix([[
#     (19 * s - 19 * t + 32) / 64,
#     (2 *s * t + 5 * s + 5 * t - 6) / 8,
# ]])
SURFACE19 = bezier.Surface(np.array([
    [0.5, -0.75],
    [0.6484375, -0.4375],
    [0.796875, -0.125],
    [0.3515625, -0.4375],
    [0.5, 0.0],
    [0.203125, -0.125],
]))
# F20 = sympy.Matrix([[
#     -(4 * s + t - 3) / 4,
#     -(16 * s**2 + 16 * s * t - 8 * s + 27 * t - 3) / 32,
# ]])
SURFACE20 = bezier.Surface(np.array([
    [0.75, 0.09375],
    [0.25, 0.21875],
    [-0.25, -0.15625],
    [0.625, -0.328125],
    [0.125, -0.453125],
    [0.5, -0.75],
]))
# F21 = sympy.Matrix([[
#     (t - s) / 2,
#     -(s**2 + s * t + 2 * s + 6 * t) / 8,
# ]])
SURFACE21 = bezier.Surface(np.array([
    [0.0, 0.0],
    [-0.25, -0.125],
    [-0.5, -0.375],
    [0.25, -0.375],
    [0.0, -0.5625],
    [0.5, -0.75],
]))
# F22 = sympy.Matrix([[(4 * s + 1) / 4, (8 * t - 1) / 8]])
SURFACE22 = bezier.Surface(np.array([
    [0.25, -0.125],
    [1.25, -0.125],
    [0.25, 0.875],
]))
# F23 = sympy.Matrix([[
#     -s**2 - s * t - 2 * t + 1, -s * (s + t - 2),
# ]])
SURFACE23 = bezier.Surface(np.array([
    [1.0, 0.0],
    [1.0, 1.0],
    [0.0, 1.0],
    [0.0, 0.0],
    [-0.5, 0.5],
    [-1.0, 0.0],
]))
# F24 = sympy.Matrix([[
#     -(3 * s**2 + 3 * s * t + 3 * t - 2) / 2,
#     -(6 * s**2 + 6 * s * t - 12 * s - t) / 4,
# ]])
SURFACE24 = bezier.Surface(np.array([
    [1.0, 0.0],
    [1.0, 1.5],
    [-0.5, 1.5],
    [0.25, 0.125],
    [-0.5, 0.875],
    [-0.5, 0.25],
]))
SURFACE25 = bezier.Surface(np.array([
    [-1.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.75],
]))
SURFACE26 = bezier.Surface(np.array([
    [-1.25, -0.25],
    [0.0, 0.25],
    [1.25, -0.25],
    [-0.375, 0.90625],
    [0.375, 0.90625],
    [0.0, 1.9375],
]))
SURFACE27 = bezier.Surface(np.array([
    [-1.25, 0.0],
    [1.25, 0.0],
    [0.0, 2.1875],
]))
SURFACE28 = bezier.Surface(np.array([
    [-1.0, 0.125],
    [0.0, -0.125],
    [1.0, 0.125],
    [-0.75, 1.1875],
    [0.75, 1.1875],
    [0.0, 1.875],
]))
# F29 = sympy.Matrix([[
#     (21 * s - 21 * t + 4) / 8,
#     (21 * s + 21 * t - 5) / 32,
# ]])
SURFACE29 = bezier.Surface(np.array([
    [0.5, -0.15625],
    [3.125, 0.5],
    [-2.125, 0.5],
]))
# F30 = sympy.Matrix([[
#     -(42 * s**2 - 10 * s * t - 33 * t**2 + 16 * s + 128 * t - 96) / 128,
#     (3 * s**2 + 24 * s * t + 11 * t**2 - 18 * t + 25) / 32,
# ]])
SURFACE30 = bezier.Surface(np.array([
    [0.75, 0.78125],
    [0.6875, 0.78125],
    [0.296875, 0.875],
    [0.25, 0.5],
    [0.2265625, 0.875],
    [0.0078125, 0.5625],
]))
# F31 = sympy.Matrix([[
#     -(s**2 + 38 * s * t + 44 * t**2 + 40 * s + 8 * t - 62) / 64,
#     (2 * s**2 + 12 * s * t + t**2 - 24 * s - 56 * t + 62) / 64,
# ]])
SURFACE31 = bezier.Surface(np.array([
    [0.96875, 0.96875],
    [0.65625, 0.78125],
    [0.328125, 0.625],
    [0.90625, 0.53125],
    [0.296875, 0.4375],
    [0.15625, 0.109375],
]))
# F32 = sympy.Matrix([[(8 - 9 * t) / 8, 2 * s + t - 1]])
SURFACE32 = bezier.Surface(np.array([
    [1.0, -1.0],
    [1.0, 1.0],
    [-0.125, 0.0],
]))
# F33 = sympy.Matrix([[(1 - 9 * s - 9 * t) / 8, s - t]])
SURFACE33 = bezier.Surface(np.array([
    [0.125, 0.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
]))
# F34 = sympy.Matrix([[s * (t + 1), t * (s + 1)]])
SURFACE34 = bezier.Surface(np.array([
    [0.0, 0.0],
    [0.5, 0.0],
    [1.0, 0.0],
    [0.0, 0.5],
    [1.0, 1.0],
    [0.0, 1.0],
]))
# F35 = sympy.Matrix([[s * t + t**2 - 2 * t + 2, s * t + t**2 + s]])
SURFACE35 = bezier.Surface(np.array([
    [2.0, 0.0],
    [2.0, 0.5],
    [2.0, 1.0],
    [1.0, 0.0],
    [1.5, 1.0],
    [1.0, 1.0],
]))
# F36 = sympy.Matrix([[
#     -(4 * s**2 + 6 * s * t - 3 * t**2 - 2 * s + 24 * t - 24) / 16,
#     (3 * s * t + 3 * t**2 + 8 * s - 3 * t) / 8,
# ]])
SURFACE36 = bezier.Surface(np.array([
    [1.5, 0.0],
    [1.5625, 0.5],
    [1.375, 1.0],
    [0.75, -0.1875],
    [0.625, 0.5],
    [0.1875, 0.0],
]))
# F37 = sympy.Matrix([[s - t, s + t]])
SURFACE37 = bezier.Surface(np.array([
    [0.0, 0.0],
    [1.0, 1.0],
    [-1.0, 1.0],
]))
# F38 = sympy.Matrix([[
#     (s - t) * (s + t + 9) / 8,
#     (s**2 - 14 * s * t + t**2 + 17 * s + 17 * t - 8) / 8,
# ]])
SURFACE38 = bezier.Surface(np.array([
    [0.0, -1.0],
    [0.5625, 0.0625],
    [1.25, 1.25],
    [-0.5625, 0.0625],
    [0.0, 0.25],
    [-1.25, 1.25],
]))


def surface_surface_check(surface1, surface2):
    if not CONFIG.running:
        return

    ax = surface1.plot(64)
    surface2.plot(64, ax=ax)
    plt.axis('scaled')
    plt.title(CONFIG.current_test)
    plt.show()


def test_surfaces1_and_3():
    surface_surface_check(SURFACE1, SURFACE3)


def test_surfaces15_and_22():
    surface_surface_check(SURFACE15, SURFACE22)


def test_surfaces1_and_2():
    surface_surface_check(SURFACE1, SURFACE2)


def test_surfaces10_and_20():
    surface_surface_check(SURFACE10, SURFACE20)


def test_surfaces10_and_21():
    surface_surface_check(SURFACE10, SURFACE21)


def test_surfaces3_and_4():
    surface_surface_check(SURFACE3, SURFACE4)


def test_surfaces1_and_29():
    surface_surface_check(SURFACE1, SURFACE29)


def test_surfaces3_and_5():
    surface_surface_check(SURFACE3, SURFACE5)


def test_surfaces15_and_16():
    surface_surface_check(SURFACE15, SURFACE16)


def test_surfaces23_and_24():
    surface_surface_check(SURFACE23, SURFACE24)


def test_surfaces25_and_26():
    surface_surface_check(SURFACE25, SURFACE26)


def test_surfaces27_and_28():
    surface_surface_check(SURFACE27, SURFACE28)


def test_surfaces6_and_7():
    surface_surface_check(SURFACE6, SURFACE7)


def test_surfaces8_and_9():
    surface_surface_check(SURFACE8, SURFACE9)


def test_surfaces4_and_10():
    surface_surface_check(SURFACE4, SURFACE10)


def test_surfaces11_and_12():
    surface_surface_check(SURFACE11, SURFACE12)


def test_surfaces3_and_13():
    surface_surface_check(SURFACE3, SURFACE13)


def test_surfaces10_and_19():
    surface_surface_check(SURFACE10, SURFACE19)


def test_surfaces3_and_14():
    surface_surface_check(SURFACE3, SURFACE14)


def test_surfaces17_and_18():
    surface_surface_check(SURFACE17, SURFACE18)


def test_surfaces30_and_31():
    surface_surface_check(SURFACE30, SURFACE31)


def test_surfaces32_and_33():
    surface_surface_check(SURFACE32, SURFACE33)


def test_surfaces34_and_35():
    surface_surface_check(SURFACE34, SURFACE35)


def test_surfaces15_and_36():
    surface_surface_check(SURFACE15, SURFACE36)


def test_surfaces37_and_38():
    surface_surface_check(SURFACE37, SURFACE38)


if __name__ == '__main__':
    CONFIG.run(globals())
