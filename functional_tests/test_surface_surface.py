# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

import collections

import matplotlib.pyplot as plt
import numpy as np
import pytest
import six

import bezier
from bezier import _plot_helpers

import runtime_utils


PARALLEL_FAILURE = ('Line segments parallel.',)
BAD_TANGENT = (
    'Curves moving in opposite direction but define '
    'overlapping arcs.')
TANGENT_FAILURE = 'The number of candidate intersections is too high.'
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
# F30Q = sympy.Matrix([[(7 * s - 2) / 8, (7 * t - 2) / 8]])
# NOTE: This is a degree elevated linear surface.
SURFACE30Q = bezier.Surface(np.array([
    [-0.25, -0.25],
    [0.1875, -0.25],
    [0.625, -0.25],
    [-0.25, 0.1875],
    [0.1875, 0.1875],
    [-0.25, 0.625],
]))
# F31Q = sympy.Matrix([[
#     -(2 * t**2 + 2 * s * t - 8 * s - t) / 8,
#     s + t - 1,
# ]])
SURFACE31Q = bezier.Surface(np.array([
    [0.0, -1.0],
    [0.5, -0.5],
    [1.0, 0.0],
    [0.0625, -0.5],
    [0.4375, 0.0],
    [-0.125, 0.0],
]))


Intersected = collections.namedtuple(
    'Intersected',
    ['start_vals', 'end_vals', 'nodes', 'edge_pairs'])


def make_plots(surface1, surface2, intersections, failed=True):
    if not CONFIG.running:
        return

    ax = surface1.plot(64)
    surface2.plot(64, ax=ax)

    color = None
    for intersection in intersections:
        intersection.plot(64, color=color, ax=ax)
        # Color is (R,G,B,A) but we just want (R,G,B).
        color = ax.patches[-1].get_facecolor()[:3]

    plt.axis('scaled')
    _plot_helpers.add_plot_boundary(ax)

    if CONFIG.save_plot:
        CONFIG.save_fig()
    else:
        if failed:
            plt.title(CONFIG.current_test + ': failed')
        else:
            plt.title(CONFIG.current_test)
        plt.show()

    plt.close(ax.figure)


def make_curved_polygon(surface1, surface2,
                        start_vals, end_vals, edge_pairs):
    base_edges = (
        surface1._get_edges(),
        surface2._get_edges(),
    )

    info = six.moves.zip(edge_pairs, start_vals, end_vals)
    edges = []
    for edge_pair, start_val, end_val in info:
        surf_index, edge_index = edge_pair
        base_edge = base_edges[surf_index][edge_index]
        edges.append(base_edge.specialize(start_val, end_val))

    return bezier.CurvedPolygon(*edges)


def surface_surface_check(surface1, surface2,
                          start_vals, end_vals, nodes, edge_pairs):
    intersected = Intersected(start_vals, end_vals, nodes, edge_pairs)
    surface_surface_check_multi(surface1, surface2, intersected)


def surface_surface_check_multi(surface1, surface2, *all_intersected):
    # pylint: disable=too-many-locals
    assert surface1.is_valid
    assert surface2.is_valid

    intersections = surface1.intersect(surface2)
    assert len(intersections) == len(all_intersected)
    edges = (
        surface1._get_edges(),
        surface2._get_edges(),
    )

    for intersection, intersected in six.moves.zip(
            intersections, all_intersected):
        assert isinstance(intersection, bezier.CurvedPolygon)

        start_vals = intersected.start_vals
        end_vals = intersected.end_vals
        nodes = intersected.nodes
        edge_pairs = intersected.edge_pairs

        info = six.moves.zip(
            intersection._edges, edge_pairs, start_vals, end_vals, nodes)
        num_edges = len(intersection._edges)
        assert num_edges == len(edge_pairs)
        assert num_edges == len(start_vals)
        assert num_edges == len(end_vals)
        assert num_edges == len(nodes)
        for edge, edge_pair, start_val, end_val, node in info:
            surf_index, edge_index = edge_pair
            expected = edges[surf_index][edge_index]
            assert expected is edge.root

            CONFIG.assert_close(edge.start, start_val)
            CONFIG.assert_close(edge.end, end_val)
            CONFIG.assert_close(edge._nodes[0, 0], node[0])
            CONFIG.assert_close(edge._nodes[0, 1], node[1])

    make_plots(surface1, surface2, intersections, failed=False)
    # pylint: enable=too-many-locals


def test_surfaces1Q_and_3Q():
    _, s_val1 = runtime_utils.real_roots([36, -48, 4, 200, -131])
    s_val2, = runtime_utils.real_roots([49, 91, 320, -244])

    t_val1, _ = runtime_utils.real_roots([9, -18, 5, -28, 12])
    t_val2, = runtime_utils.real_roots([49, 63, 88, -128])
    start_vals = np.array([s_val2, 0.0, 0.0, t_val1])
    end_vals = np.array([1.0, 1.0, s_val1, t_val2])

    x_val1 = 0.5 * (1.0 - s_val1) * (s_val1 + 2.0)
    y_val1 = 0.5 * s_val1 * (3.0 - s_val1)
    x_val2 = 0.5 * s_val2 * (1.0 - s_val2)
    y_val2 = 1.0 - s_val2
    nodes = np.array([
        [x_val2, y_val2],
        [0.0, 0.0],
        [1.0, 0.0],
        [x_val1, y_val1],
    ])
    edge_pairs = (
        (0, 2),
        (1, 0),
        (0, 1),
        (1, 2),
    )
    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(45):
        surface_surface_check(SURFACE1Q, SURFACE3Q,
                              start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_3L():
    start_vals = np.array([0.125, 0.25, 0.0])
    end_vals = np.array([0.875, 1.0, 0.75])

    nodes = np.array([
        [0.25, 0.75],
        [0.25, 0.0],
        [1.0, 0.0],
    ])
    edge_pairs = (
        (1, 2),
        (0, 0),
        (0, 1),
    )
    surface_surface_check(SURFACE1L, SURFACE3L,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1Q_and_2Q():
    # pylint: disable=too-many-locals
    s_val1, _ = runtime_utils.real_roots([3, -21, 10])
    _, s_val2 = runtime_utils.real_roots([12, -24, 0, 140, -105])
    s_val3, _ = runtime_utils.real_roots([12, -72, 56, -100, 23])
    _, s_val4 = runtime_utils.real_roots([12, 24, -88, 156, -81])
    _, s_val5 = runtime_utils.real_roots([9, -54, 213, -12, -96])
    s_val6, _ = runtime_utils.real_roots([12, -24, 24, -140, 23])

    _, t_val1 = runtime_utils.real_roots([9, 39, -38])
    t_val2, _ = runtime_utils.real_roots([9, -18, -3, -116, 20])
    _, t_val3 = runtime_utils.real_roots([9, 30, -119, 272, -128])
    t_val4, _ = runtime_utils.real_roots([9, -66, 25, -160, 64])
    _, t_val5 = runtime_utils.real_roots([9, -66, 181, 36, -44])
    t_val6, _ = runtime_utils.real_roots([9, -18, -3, -116, 84])
    start_vals = np.array(
        [s_val6, t_val5, s_val1, t_val2, s_val3, t_val4])
    end_vals = np.array(
        [s_val5, t_val1, s_val2, t_val3, s_val4, t_val6])

    x_val1 = 0.5 * s_val6 * (1.0 - s_val6)
    x_val2 = 0.5 * s_val5 * (1.0 - s_val5)
    x_val5 = 0.5 * (1.0 - s_val3) * (s_val3 + 2.0)
    x_val6 = 0.5 * (1.0 - s_val4) * (s_val4 + 2.0)

    y_val1 = 1.0 - s_val6
    y_val2 = 1.0 - s_val5
    y_val3 = 0.5 * s_val1 * (s_val1 - 1.0)
    y_val4 = 0.5 * s_val2 * (s_val2 - 1.0)
    y_val5 = 0.5 * s_val3 * (3.0 - s_val3)
    y_val6 = 0.5 * s_val4 * (3.0 - s_val4)

    nodes = np.array([
        [x_val1, y_val1],
        [x_val2, y_val2],
        [s_val1, y_val3],
        [s_val2, y_val4],
        [x_val5, y_val5],
        [x_val6, y_val6],
    ])
    edge_pairs = (
        (0, 2),
        (1, 0),
        (0, 0),
        (1, 1),
        (0, 1),
        (1, 2),
    )
    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(9):
        surface_surface_check(SURFACE1Q, SURFACE2Q,
                              start_vals, end_vals, nodes, edge_pairs)
    # pylint: enable=too-many-locals


def test_surfaces10Q_and_18Q():
    start_vals = np.array([0.0, 0.25, 0.0])
    end_vals = np.array([1.0, 1.0, 1.0])

    nodes = np.array([
        [0.5, -0.75],
        [0.75, 0.09375],
        [0.0, 0.0],
    ])
    edge_pairs = (
        (1, 2),
        (0, 0),
        (0, 1),
    )
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check(SURFACE10Q, SURFACE18Q,
                              start_vals, end_vals, nodes, edge_pairs)

    assert str(exc_info.value).startswith(TANGENT_FAILURE)
    intersection = make_curved_polygon(
        SURFACE10Q, SURFACE18Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE10Q, SURFACE18Q, [intersection])


def test_surfaces10Q_and_19Q():
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check_multi(SURFACE10Q, SURFACE19Q)

    assert exc_info.value.args == PARALLEL_FAILURE
    make_plots(SURFACE10Q, SURFACE19Q, [])


def test_surfaces3Q_and_4Q():
    _, s_val2 = runtime_utils.real_roots([25, -130, 321, -88, -96])
    _, s_val3 = runtime_utils.real_roots([49, 14, 145, 1008, -468])

    t_val2, _ = runtime_utils.real_roots([100, 360, 712, -2988, 169])
    _, t_val3 = runtime_utils.real_roots([49, -532, 412, 37200, -26352])
    start_vals = np.array([s_val3, t_val2, 0.0])
    end_vals = np.array([s_val2, 1.0, t_val3])

    x_val2 = 0.125 * (s_val2 - 1.0) * (5.0 * s_val2 - 8.0)
    x_val3 = 0.125 * (s_val3 - 1.0) * (5.0 * s_val3 - 8.0)
    y_val2 = 0.125 * (1.0 - s_val2) * (7.0 * s_val2 + 8.0)
    y_val3 = 0.125 * (1.0 - s_val3) * (7.0 * s_val3 + 8.0)
    nodes = np.array([
        [x_val3, y_val3],
        [x_val2, y_val2],
        [1.0, 0.25],
    ])
    edge_pairs = (
        (0, 2),
        (1, 0),
        (1, 1),
    )
    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(32):
        surface_surface_check(SURFACE3Q, SURFACE4Q,
                              start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1Q_and_5L():
    s_val4, _ = runtime_utils.real_roots([1, -3, 1])
    t_val4, _ = runtime_utils.real_roots([1764, -3108, 1049])
    start_vals = np.array([0.5, 0.0, 0.0, t_val4])
    end_vals = np.array([1.0, 1.0, s_val4, 4.0 / 7.0])

    x_val4 = 0.5 * (1.0 - s_val4) * (s_val4 + 2.0)
    nodes = np.array([
        [0.125, 0.5],
        [0.0, 0.0],
        [1.0, 0.0],
        [x_val4, 0.5],
    ])
    edge_pairs = (
        (0, 2),
        (0, 0),
        (0, 1),
        (1, 1),
    )
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check(SURFACE1Q, SURFACE5L,
                              start_vals, end_vals, nodes, edge_pairs)

    assert str(exc_info.value).startswith(TANGENT_FAILURE)
    intersection = make_curved_polygon(
        SURFACE1Q, SURFACE5L,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE1Q, SURFACE5L, [intersection])


def test_surfaces3Q_and_5Q():
    s_val3, _ = runtime_utils.real_roots([25, -130, 167, -302, 57])
    s_val4, _ = runtime_utils.real_roots([25, -130, 901, -1212, 279])

    _, t_val3 = runtime_utils.real_roots([25, -20, -1064, 7800, -6012])
    t_val4, _ = runtime_utils.real_roots([25, -2340, 58908, -105840, 11664])
    start_vals = np.array([t_val4, 0.0, 0.0, s_val3])
    end_vals = np.array([1.0, 1.0, t_val3, s_val4])

    x_val3 = 0.125 * (s_val3 - 1.0) * (5.0 * s_val3 - 8.0)
    x_val4 = 0.125 * (s_val4 - 1.0) * (5.0 * s_val4 - 8.0)
    y_val3 = 0.125 * (1.0 - s_val3) * (7.0 * s_val3 + 8.0)
    y_val4 = 0.125 * (1.0 - s_val4) * (7.0 * s_val4 + 8.0)
    nodes = np.array([
        [x_val4, y_val4],
        [0.25, 0.09375],
        [1.125, 0.375],
        [x_val3, y_val3],
    ])
    edge_pairs = (
        (1, 2),
        (1, 0),
        (1, 1),
        (0, 2),
    )
    surface_surface_check(SURFACE3Q, SURFACE5Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_2L():
    start_vals = np.array([7.59375, 1.0, 3.0, 4.5, 0.0]) / 9.0
    end_vals = np.array([27.0, 19.0, 11.8125, 27.0, 13.5]) / 27.0

    nodes = np.array([
        [0.0, 0.46875],
        [0.0, 0.0],
        [2.0, 1.0],
        [1.6875, 1.3125],
        [0.375, 1.125],
    ]) / 3.0
    edge_pairs = (
        (0, 2),
        (1, 0),
        (0, 1),
        (1, 1),
        (1, 2),
    )
    surface_surface_check(SURFACE1L, SURFACE2L,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces20Q_and_21Q():
    start_vals = np.array([0.0, 0.0, 4.0 / 5.0, 0.0])
    end_vals = np.array([1.0, 0.5, 1.0, 1.0])

    nodes = np.array([
        [1.0, 0.0],
        [0.0, 1.0],
        [-0.5, 0.5],
        [-0.5, 0.25],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 2),
    )
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check(SURFACE20Q, SURFACE21Q,
                              start_vals, end_vals, nodes, edge_pairs)

    assert exc_info.value.args == PARALLEL_FAILURE
    intersection = make_curved_polygon(
        SURFACE20Q, SURFACE21Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE20Q, SURFACE21Q, [intersection])


def test_surfaces4L_and_22Q():
    start_vals = np.array([0.0, 0.0, 0.0])
    end_vals = np.array([1.0, 1.0, 1.0])

    nodes = np.array([
        [-1.25, 0.0],
        [1.25, 0.0],
        [0.0, 2.1875],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (0, 2),
    )
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check(SURFACE4L, SURFACE22Q,
                              start_vals, end_vals, nodes, edge_pairs)

    assert str(exc_info.value).startswith(TANGENT_FAILURE)
    intersection = make_curved_polygon(
        SURFACE4L, SURFACE22Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE4L, SURFACE22Q, [intersection])


def test_surfaces4L_and_23Q():
    start_vals = np.array([0.0, 0.0, 0.0])
    end_vals = np.array([1.0, 1.0, 1.0])

    nodes = np.array([
        [-1.0, 0.125],
        [1.0, 0.125],
        [0.0, 1.875],
    ])
    edge_pairs = (
        (1, 0),
        (1, 1),
        (1, 2),
    )
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check(SURFACE4L, SURFACE23Q,
                              start_vals, end_vals, nodes, edge_pairs)

    assert str(exc_info.value).startswith(TANGENT_FAILURE)
    intersection = make_curved_polygon(
        SURFACE4L, SURFACE23Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE4L, SURFACE23Q, [intersection])


def test_surfaces6Q_and_7Q():
    # pylint: disable=too-many-locals
    s_val3, _ = runtime_utils.real_roots([1, -13, 2])
    _, s_val4 = runtime_utils.real_roots([7, 5, -10])
    s_val5, s_val6 = runtime_utils.real_roots([4, 120, 1592, -1908, 489])
    s_val7, _ = runtime_utils.real_roots([64, -1232, 297])
    _, s_val8 = runtime_utils.real_roots([576, 784, -871])

    t_val3, _ = runtime_utils.real_roots([1, -102, 25])
    _, t_val4 = runtime_utils.real_roots([49, 68, -76])
    t_val6, t_val5 = runtime_utils.real_roots([4, -104, 1504, -1548, 369])
    t_val7, _ = runtime_utils.real_roots([4, -20, 3])
    _, t_val8 = runtime_utils.real_roots([12, 4, -13])

    start_vals1 = np.array([s_val8, 3.0 / 17.0, s_val3, t_val5])
    end_vals1 = np.array([14.0 / 17.0, t_val3, s_val5, t_val8])
    start_vals2 = np.array([t_val7, s_val6, t_val4, 3.0 / 17.0])
    end_vals2 = np.array([t_val6, s_val4, 14.0 / 17.0, s_val7])

    x_val3 = 0.25 * s_val3 * (3.0 + s_val3)
    x_val4 = 0.25 * s_val4 * (3.0 + s_val4)
    x_val5 = 0.25 * s_val5 * (3.0 + s_val5)
    x_val6 = 0.25 * s_val6 * (3.0 + s_val6)
    y_val3 = 0.75 * s_val3 * (1.0 - s_val3)
    y_val4 = 0.75 * s_val4 * (1.0 - s_val4)
    y_val5 = 0.75 * s_val5 * (1.0 - s_val5)
    y_val6 = 0.75 * s_val6 * (1.0 - s_val6)

    nodes1 = np.array([
        [0.5 - 0.5 * s_val8, 1.0 - s_val8],
        [3.0 / 34.0, 3.0 / 17.0],
        [x_val3, y_val3],
        [x_val5, y_val5],
    ])
    edge_pairs1 = (
        (0, 2),
        (1, 1),
        (0, 0),
        (1, 0),
    )

    nodes2 = np.array([
        [1.0 - 0.5 * s_val7, s_val7],
        [x_val6, y_val6],
        [x_val4, y_val4],
        [31.0 / 34.0, 3.0 / 17.0],
    ])
    edge_pairs2 = (
        (1, 0),
        (0, 0),
        (1, 2),
        (0, 1),
    )
    intersected1 = Intersected(start_vals1, end_vals1, nodes1, edge_pairs1)
    intersected2 = Intersected(start_vals2, end_vals2, nodes2, edge_pairs2)
    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(25):
        surface_surface_check_multi(SURFACE6Q, SURFACE7Q,
                                    intersected1, intersected2)
    # pylint: enable=too-many-locals


def test_surfaces8Q_and_9Q():
    s_val2, s_val3 = runtime_utils.real_roots([28, -24, 1])
    t_val3, t_val2 = runtime_utils.real_roots([28, -32, 5])
    start_vals = np.array([11.0 / 14.0, s_val2, 1.0 / 14.0, s_val3])
    end_vals = np.array([t_val2, 3.0 / 14.0, t_val3, 13.0 / 14.0])

    x_val2 = (1.0 - s_val2) * (1.0 + 7.0 * s_val2)
    x_val3 = (1.0 - s_val3) * (1.0 + 7.0 * s_val3)
    y_val2 = s_val2 * (8.0 - 7.0 * s_val2)
    y_val3 = s_val3 * (8.0 - 7.0 * s_val3)
    nodes = np.array([
        [15.0 / 28.0, 39.0 / 28.0],
        [x_val2, y_val2],
        [55.0 / 28.0, 39.0 / 28.0],
        [x_val3, y_val3],
    ])
    edge_pairs = (
        (1, 2),
        (0, 1),
        (1, 2),
        (0, 1),
    )
    surface_surface_check(SURFACE8Q, SURFACE9Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces4Q_and_10Q():
    surface_surface_check_multi(SURFACE4Q, SURFACE10Q)


def test_surfaces11Q_and_12Q():
    s_val1, s_val2 = runtime_utils.real_roots([8, -8, 1])
    start_vals = np.array([s_val1, s_val1])
    end_vals = np.array([s_val2, s_val2])

    nodes = np.array([
        [s_val2, 0.125],
        [s_val1, 0.125],
    ])
    edge_pairs = (
        (1, 0),
        (0, 0),
    )
    surface_surface_check(SURFACE11Q, SURFACE12Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces3Q_and_13Q():
    start_vals = np.array([0.0, 0.0, 0.0])
    end_vals = np.array([1.0, 1.0, 1.0])

    nodes = np.array([
        [0.25, 0.15625],
        [0.75, 0.15625],
        [0.5, 0.625],
    ])
    edge_pairs = (
        (1, 0),
        (1, 1),
        (1, 2),
    )
    surface_surface_check(SURFACE3Q, SURFACE13Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces10Q_and_17Q():
    start_vals = np.array([0.0, 0.0, 0.0])
    end_vals = np.array([1.0, 1.0, 1.0])

    nodes = np.array([
        [0.5, -0.75],
        [0.796875, -0.125],
        [0.203125, -0.125],
    ])
    edge_pairs = (
        (1, 0),
        (1, 1),
        (1, 2),
    )
    surface_surface_check(SURFACE10Q, SURFACE17Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces3Q_and_14Q():
    start_vals = np.array([0.0, 0.0, 0.0])
    end_vals = np.array([1.0, 1.0, 1.0])

    nodes = np.array([
        [0.25, 0.21875],
        [0.75, 0.21875],
        [0.5, 0.6875],
    ])
    edge_pairs = (
        (1, 0),
        (1, 1),
        (1, 2),
    )

    surface_surface_check(SURFACE3Q, SURFACE14Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces15Q_and_16Q():
    # pylint: disable=too-many-locals
    s_val4, _ = runtime_utils.real_roots([49, -120, 32])
    _, s_val5 = runtime_utils.real_roots([1, 70, -39])
    s_val6, _ = runtime_utils.real_roots([2, -18, 1])
    _, s_val7 = runtime_utils.real_roots([14, 2, -9])

    t_val4, _ = runtime_utils.real_roots([14, -30, 7])
    _, t_val5 = runtime_utils.real_roots([2, 14, -15])
    t_val6, _ = runtime_utils.real_roots([1, -72, 32])
    _, t_val7 = runtime_utils.real_roots([49, 22, -39])

    start_vals1 = np.array([0.25, t_val4, 0.5, t_val7])
    end_vals1 = np.array([s_val4, 0.5, s_val7, 0.75])
    start_vals2 = np.array([s_val6, 0.5, s_val5, 0.25])
    end_vals2 = np.array([0.5, t_val5, 0.75, t_val6])

    nodes1 = np.array([
        [0.5, 0.4375],
        [2.0 * s_val4, 1.75 * s_val4],
        [0.5, 1.0],
        [2.0 * s_val7 * (1.0 - s_val7), 2.0 - 2.0 * s_val7],
    ])
    edge_pairs1 = (
        (0, 0),
        (1, 0),
        (0, 2),
        (1, 2),
    )

    nodes2 = np.array([
        [2.0 * s_val6 * (1.0 - s_val6), 2.0 - 2.0 * s_val6],
        [0.5, 1.0],
        [2.0 - 2.0 * s_val5, 1.75 + 0.25 * s_val5],
        [0.5, 1.9375],
    ])
    edge_pairs2 = (
        (0, 2),
        (1, 0),
        (0, 1),
        (1, 1),
    )
    intersected1 = Intersected(start_vals1, end_vals1, nodes1, edge_pairs1)
    intersected2 = Intersected(start_vals2, end_vals2, nodes2, edge_pairs2)

    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check_multi(SURFACE15Q, SURFACE16Q,
                                    intersected1, intersected2)

    assert exc_info.value.args == (BAD_TANGENT,)
    intersection1 = make_curved_polygon(
        SURFACE15Q, SURFACE16Q,
        start_vals1, end_vals1, edge_pairs1)
    intersection2 = make_curved_polygon(
        SURFACE15Q, SURFACE16Q,
        start_vals2, end_vals2, edge_pairs2)
    make_plots(SURFACE15Q, SURFACE16Q, [intersection1, intersection2])
    # pylint: enable=too-many-locals


def test_surfaces24Q_and_25Q():
    _, _, s_val1, _ = runtime_utils.real_roots(
        [81, 72, -4640, -1168, 1036])
    _, s_val2 = runtime_utils.real_roots(
        [121, -14740, 618410, -9692, -153359])
    _, t_val1, _, _ = runtime_utils.real_roots(
        [27, -1116, 12020, -10224, 2256])
    _, t_val2 = runtime_utils.real_roots(
        [11, -1232, 132116, 315936, -31348])
    start_vals = np.array([s_val2, 0.0, t_val1, 0.0])
    end_vals = np.array([1.0, s_val1, 1.0, t_val2])

    x_val1 = 0.015625 * (4.0 - 3.0 * s_val1) * (7.0 * s_val1 + 12.0)
    y_val1 = 0.03125 * (3.0 * s_val1 * s_val1 + 25.0)
    x_val2 = 0.0078125 * (33.0 * s_val2 * s_val2 + 62.0 * s_val2 + 1.0)
    y_val2 = 0.03125 * (11.0 * s_val2 * s_val2 - 4.0 * s_val2 + 18.0)
    nodes = np.array([
        [x_val2, y_val2],
        [0.75, 0.78125],
        [x_val1, y_val1],
        [0.328125, 0.625],
    ])
    edge_pairs = (
        (0, 2),
        (0, 0),
        (1, 0),
        (1, 1),
    )

    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(22):
        surface_surface_check(SURFACE24Q, SURFACE25Q,
                              start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_6L():
    start_vals = np.array([0.75, 0.0, 0.75, 0.0])
    end_vals = np.array([1.0, 0.25, 1.0, 0.25])

    nodes = np.array([
        [0.0, 0.25],
        [0.0, 0.0],
        [0.25, 0.0],
        [0.25, 0.25],
    ])
    edge_pairs = (
        (0, 2),
        (0, 0),
        (1, 2),
        (1, 0),
    )
    surface_surface_check(SURFACE1L, SURFACE6L,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces26Q_and_27Q():
    surface_surface_check_multi(SURFACE26Q, SURFACE27Q)


def test_surfaces1L_and_28Q():
    _, s_val3 = runtime_utils.real_roots([5, 30, -13])
    t_val3, _ = runtime_utils.real_roots([5, -40, 22])
    start_vals = np.array([t_val3, 0.1875, 0.0])
    end_vals = np.array([1.0, 1.0, s_val3])

    nodes = np.array([
        [1.0 - s_val3, s_val3],
        [0.1875, 0.0],
        [1.0, 0.0],
    ])
    edge_pairs = (
        (1, 1),
        (0, 0),
        (0, 1),
    )
    surface_surface_check(SURFACE1L, SURFACE28Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_29Q():
    s_val1, s_val2 = runtime_utils.real_roots([128, -128, 7])
    t_val1, t_val2 = runtime_utils.real_roots([8, -8, 1])
    start_vals = np.array([s_val2, 0.0, 0.0, 0.0, t_val1])
    end_vals = np.array([1.0, 1.0, 1.0, s_val1, t_val2])

    nodes = np.array([
        [1.0 - s_val2, s_val2],
        [0.0, 1.0],
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0 - s_val1, s_val1],
    ])
    edge_pairs = (
        (0, 1),
        (0, 2),
        (0, 0),
        (0, 1),
        (1, 1),
    )
    surface_surface_check(SURFACE1L, SURFACE29Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces30Q_and_31Q():
    start_vals = np.array([13.0 / 56.0, 0.0, 5.0 / 9.0, 0.0])
    end_vals = np.array([1.0, 2.0 / 7.0, 1.0, 0.25])

    nodes = np.array([
        [-0.046875, -0.25],
        [0.625, -0.25],
        [0.375, 0.0],
        [-0.125, 0.0],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 2),
    )

    surface_surface_check(SURFACE30Q, SURFACE31Q,
                          start_vals, end_vals, nodes, edge_pairs)


if __name__ == '__main__':
    CONFIG.run(globals())
