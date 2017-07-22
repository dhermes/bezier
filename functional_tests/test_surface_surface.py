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

from __future__ import absolute_import
import collections
import contextlib
import operator

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
import numpy as np
import pytest
import six

import bezier
from bezier import _implicitization
from bezier import _plot_helpers
from bezier import curve

import runtime_utils


ALGEBRAIC = curve.IntersectionStrategy.algebraic
GEOMETRIC = curve.IntersectionStrategy.geometric
SURFACES, INTERSECTIONS = runtime_utils.surface_intersections_info()
STRATEGY = GEOMETRIC
PARALLEL_FAILURE = ('Line segments parallel.',)
BAD_TANGENT = (
    'Curves moving in opposite direction but define '
    'overlapping arcs.')
BAD_TANGENT = (BAD_TANGENT,)
TANGENT_FAILURE = 'The number of candidate intersections is too high.'
WIGGLES = {
    1: 46,
    13: 19,
    32: 1013,
    33: 1013,
}
FAILED_CASES_TANGENT = {
    10: {'parallel': True},
    21: {'bad_tangent': True},
}
CONFIG = runtime_utils.Config()

SURFACE1L = SURFACES['1L'].surface
SURFACE2L = SURFACES['2L'].surface
SURFACE3L = SURFACES['3L'].surface
SURFACE4L = SURFACES['4L'].surface
SURFACE5L = SURFACES['5L'].surface
SURFACE6L = SURFACES['6L'].surface
SURFACE7L = SURFACES['7L'].surface
SURFACE8L = SURFACES['8L'].surface
SURFACE9L = SURFACES['9L'].surface
SURFACE10L = SURFACES['10L'].surface
SURFACE1Q = SURFACES['1Q'].surface
SURFACE2Q = SURFACES['2Q'].surface
SURFACE3Q = SURFACES['3Q'].surface
SURFACE4Q = SURFACES['4Q'].surface
SURFACE5Q = SURFACES['5Q'].surface
SURFACE10Q = SURFACES['10Q'].surface
SURFACE13Q = SURFACES['13Q'].surface
SURFACE14Q = SURFACES['14Q'].surface
SURFACE17Q = SURFACES['17Q'].surface
SURFACE18Q = SURFACES['18Q'].surface
SURFACE19Q = SURFACES['19Q'].surface
SURFACE22Q = SURFACES['22Q'].surface
SURFACE23Q = SURFACES['23Q'].surface
SURFACE24Q = SURFACES['24Q'].surface
SURFACE25Q = SURFACES['25Q'].surface
SURFACE26Q = SURFACES['26Q'].surface
SURFACE27Q = SURFACES['27Q'].surface
SURFACE28Q = SURFACES['28Q'].surface
SURFACE29Q = SURFACES['29Q'].surface
SURFACE30Q = SURFACES['30Q'].surface
SURFACE31Q = SURFACES['31Q'].surface


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


def curved_polygon_edges(intersection, edges):
    edges1, edges2 = edges
    all_edges = edges1 + edges2
    # Re-sort the edges to be in the same order independent of strategy.
    edge_list = intersection._edges
    edge_info = [(all_edges.index(edge.root), edge.start, edge.end)
                 for edge in edge_list]
    index = edge_info.index(min(edge_info))
    return edge_list[index:] + edge_list[:index]


def check_tangent(caught_exc, parallel=False, bad_tangent=False):
    exc_args = caught_exc.args
    if STRATEGY is GEOMETRIC:
        if parallel:
            assert exc_args == PARALLEL_FAILURE
        elif bad_tangent:
            assert exc_args == BAD_TANGENT
        else:
            assert len(exc_args) == 1
            assert exc_args[0].startswith(TANGENT_FAILURE)
    else:
        assert len(exc_args) == 2
        assert exc_args[0] == _implicitization._NON_SIMPLE_ERR


@contextlib.contextmanager
def check_tangent_manager(**kwargs):
    caught_exc = None
    try:
        yield
    except NotImplementedError as exc:
        caught_exc = exc

    assert caught_exc is not None
    check_tangent(caught_exc, **kwargs)


def check_coincident(exc_info, parallel=False):
    if STRATEGY is GEOMETRIC:
        if parallel:
            assert exc_info.value.args == PARALLEL_FAILURE
        else:
            assert str(exc_info.value).startswith(TANGENT_FAILURE)
    else:
        assert exc_info.value.args == (_implicitization._COINCIDENT_ERR,)


def surface_surface_check_multi(surface1, surface2, *all_intersected):
    # pylint: disable=too-many-locals
    assert surface1.is_valid
    assert surface2.is_valid

    intersections = surface1.intersect(surface2, strategy=STRATEGY)
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

        int_edges = curved_polygon_edges(intersection, edges)
        info = six.moves.zip(
            int_edges, edge_pairs, start_vals, end_vals, nodes)
        num_edges = len(int_edges)
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


def test_surfaces1L_and_3L():
    start_vals = np.asfortranarray([0.25, 0.0, 0.125])
    end_vals = np.asfortranarray([1.0, 0.75, 0.875])

    nodes = np.asfortranarray([
        [0.25, 0.0],
        [1.0, 0.0],
        [0.25, 0.75],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (1, 2),
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
    start_vals = np.asfortranarray(
        [s_val1, t_val2, s_val3, t_val4, s_val6, t_val5])
    end_vals = np.asfortranarray(
        [s_val2, t_val3, s_val4, t_val6, s_val5, t_val1])

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

    nodes = np.asfortranarray([
        [s_val1, y_val3],
        [s_val2, y_val4],
        [x_val5, y_val5],
        [x_val6, y_val6],
        [x_val1, y_val1],
        [x_val2, y_val2],
    ])
    edge_pairs = (
        (0, 0),
        (1, 1),
        (0, 1),
        (1, 2),
        (0, 2),
        (1, 0),
    )
    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(19):
        surface_surface_check(SURFACE1Q, SURFACE2Q,
                              start_vals, end_vals, nodes, edge_pairs)
    # pylint: enable=too-many-locals


def test_surfaces10Q_and_18Q():
    start_vals = np.asfortranarray([0.0, 0.25, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
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

    check_coincident(exc_info)
    intersection = make_curved_polygon(
        SURFACE10Q, SURFACE18Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE10Q, SURFACE18Q, [intersection])


def test_surfaces10Q_and_19Q():
    with pytest.raises(NotImplementedError) as exc_info:
        surface_surface_check_multi(SURFACE10Q, SURFACE19Q)

    check_coincident(exc_info, parallel=True)
    make_plots(SURFACE10Q, SURFACE19Q, [])


def test_surfaces3Q_and_4Q():
    _, s_val2 = runtime_utils.real_roots([25, -130, 321, -88, -96])
    _, s_val3 = runtime_utils.real_roots([49, 14, 145, 1008, -468])

    t_val2, _ = runtime_utils.real_roots([100, 360, 712, -2988, 169])
    _, t_val3 = runtime_utils.real_roots([49, -532, 412, 37200, -26352])
    start_vals = np.asfortranarray([s_val3, t_val2, 0.0])
    end_vals = np.asfortranarray([s_val2, 1.0, t_val3])

    x_val2 = 0.125 * (s_val2 - 1.0) * (5.0 * s_val2 - 8.0)
    x_val3 = 0.125 * (s_val3 - 1.0) * (5.0 * s_val3 - 8.0)
    y_val2 = 0.125 * (1.0 - s_val2) * (7.0 * s_val2 + 8.0)
    y_val3 = 0.125 * (1.0 - s_val3) * (7.0 * s_val3 + 8.0)
    nodes = np.asfortranarray([
        [x_val3, y_val3],
        [x_val2, y_val2],
        [1.0, 0.25],
    ])
    edge_pairs = (
        (0, 2),
        (1, 0),
        (1, 1),
    )
    if STRATEGY is GEOMETRIC:
        # NOTE: We require a bit more wiggle room for these roots.
        with CONFIG.wiggle(32):
            surface_surface_check(SURFACE3Q, SURFACE4Q,
                                  start_vals, end_vals, nodes, edge_pairs)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            surface_surface_check(SURFACE3Q, SURFACE4Q,
                                  start_vals, end_vals, nodes, edge_pairs)

        check_tangent(exc_info.value)


def test_surfaces1Q_and_5L():
    s_val4, _ = runtime_utils.real_roots([1, -3, 1])
    t_val4, _ = runtime_utils.real_roots([1764, -3108, 1049])
    start_vals = np.asfortranarray([0.5, 0.0, 0.0, t_val4])
    end_vals = np.asfortranarray([1.0, 1.0, s_val4, 4.0 / 7.0])

    x_val4 = 0.5 * (1.0 - s_val4) * (s_val4 + 2.0)
    nodes = np.asfortranarray([
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

    check_tangent(exc_info.value)
    intersection = make_curved_polygon(
        SURFACE1Q, SURFACE5L,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE1Q, SURFACE5L, [intersection])


def test_surfaces3Q_and_5Q():
    s_val3, _ = runtime_utils.real_roots([25, -130, 167, -302, 57])
    s_val4, _ = runtime_utils.real_roots([25, -130, 901, -1212, 279])

    _, t_val3 = runtime_utils.real_roots([25, -20, -1064, 7800, -6012])
    t_val4, _ = runtime_utils.real_roots([25, -2340, 58908, -105840, 11664])
    start_vals = np.asfortranarray([s_val3, t_val4, 0.0, 0.0])
    end_vals = np.asfortranarray([s_val4, 1.0, 1.0, t_val3])

    x_val3 = 0.125 * (s_val3 - 1.0) * (5.0 * s_val3 - 8.0)
    x_val4 = 0.125 * (s_val4 - 1.0) * (5.0 * s_val4 - 8.0)
    y_val3 = 0.125 * (1.0 - s_val3) * (7.0 * s_val3 + 8.0)
    y_val4 = 0.125 * (1.0 - s_val4) * (7.0 * s_val4 + 8.0)
    nodes = np.asfortranarray([
        [x_val3, y_val3],
        [x_val4, y_val4],
        [0.25, 0.09375],
        [1.125, 0.375],
    ])
    edge_pairs = (
        (0, 2),
        (1, 2),
        (1, 0),
        (1, 1),
    )
    surface_surface_check(SURFACE3Q, SURFACE5Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_2L():
    start_vals = np.asfortranarray([3.0, 4.5, 0.0, 7.59375, 1.0]) / 9.0
    end_vals = np.asfortranarray([11.8125, 27.0, 13.5, 27.0, 19.0]) / 27.0

    nodes = np.asfortranarray([
        [2.0, 1.0],
        [1.6875, 1.3125],
        [0.375, 1.125],
        [0.0, 0.46875],
        [0.0, 0.0],
    ]) / 3.0
    edge_pairs = (
        (0, 1),
        (1, 1),
        (1, 2),
        (0, 2),
        (1, 0),
    )
    surface_surface_check(SURFACE1L, SURFACE2L,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces4L_and_22Q():
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
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

    check_tangent(exc_info.value)
    intersection = make_curved_polygon(
        SURFACE4L, SURFACE22Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE4L, SURFACE22Q, [intersection])


def test_surfaces4L_and_23Q():
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
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

    check_tangent(exc_info.value)
    intersection = make_curved_polygon(
        SURFACE4L, SURFACE23Q,
        start_vals, end_vals, edge_pairs)
    make_plots(SURFACE4L, SURFACE23Q, [intersection])


def test_surfaces4Q_and_10Q():
    if STRATEGY is GEOMETRIC:
        surface_surface_check_multi(SURFACE4Q, SURFACE10Q)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            surface_surface_check_multi(SURFACE4Q, SURFACE10Q)

        check_tangent(exc_info.value)


def test_surfaces3Q_and_13Q():
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
        [0.25, 0.15625],
        [0.75, 0.15625],
        [0.5, 0.625],
    ])
    edge_pairs = (
        (1, 0),
        (1, 1),
        (1, 2),
    )
    if STRATEGY is GEOMETRIC:
        surface_surface_check(SURFACE3Q, SURFACE13Q,
                              start_vals, end_vals, nodes, edge_pairs)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            surface_surface_check(SURFACE3Q, SURFACE13Q,
                                  start_vals, end_vals, nodes, edge_pairs)

        check_tangent(exc_info.value)


def test_surfaces10Q_and_17Q():
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
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


def test_surfaces17Q_and_10Q():
    # NOTE: This is identical to test_surfaces10Q_and_17Q, but
    #       now the "first" surface is the one that is fully
    #       interior at the double corner.
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
        [0.5, -0.75],
        [0.796875, -0.125],
        [0.203125, -0.125],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (0, 2),
    )
    surface_surface_check(SURFACE17Q, SURFACE10Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces3Q_and_14Q():
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
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


def test_surfaces24Q_and_25Q():
    _, _, s_val1, _ = runtime_utils.real_roots(
        [81, 72, -4640, -1168, 1036])
    _, s_val2 = runtime_utils.real_roots(
        [121, -14740, 618410, -9692, -153359])
    _, t_val1, _, _ = runtime_utils.real_roots(
        [27, -1116, 12020, -10224, 2256])
    _, t_val2 = runtime_utils.real_roots(
        [11, -1232, 132116, 315936, -31348])
    start_vals = np.asfortranarray([0.0, t_val1, 0.0, s_val2])
    end_vals = np.asfortranarray([s_val1, 1.0, t_val2, 1.0])

    x_val1 = 0.015625 * (4.0 - 3.0 * s_val1) * (7.0 * s_val1 + 12.0)
    y_val1 = 0.03125 * (3.0 * s_val1 * s_val1 + 25.0)
    x_val2 = 0.0078125 * (33.0 * s_val2 * s_val2 + 62.0 * s_val2 + 1.0)
    y_val2 = 0.03125 * (11.0 * s_val2 * s_val2 - 4.0 * s_val2 + 18.0)
    nodes = np.asfortranarray([
        [0.75, 0.78125],
        [x_val1, y_val1],
        [0.328125, 0.625],
        [x_val2, y_val2],
    ])
    edge_pairs = (
        (0, 0),
        (1, 0),
        (1, 1),
        (0, 2),
    )

    # NOTE: We require a bit more wiggle room for these roots.
    with CONFIG.wiggle(35):
        surface_surface_check(SURFACE24Q, SURFACE25Q,
                              start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_6L():
    start_vals = np.asfortranarray([0.0, 0.75, 0.0, 0.75])
    end_vals = np.asfortranarray([0.25, 1.0, 0.25, 1.0])

    nodes = np.asfortranarray([
        [0.0, 0.0],
        [0.25, 0.0],
        [0.25, 0.25],
        [0.0, 0.25],
    ])
    edge_pairs = (
        (0, 0),
        (1, 2),
        (1, 0),
        (0, 2),
    )
    surface_surface_check(SURFACE1L, SURFACE6L,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces26Q_and_27Q():
    surface_surface_check_multi(SURFACE26Q, SURFACE27Q)


def test_surfaces1L_and_28Q():
    _, s_val3 = runtime_utils.real_roots([5, 30, -13])
    t_val3, _ = runtime_utils.real_roots([5, -40, 22])
    start_vals = np.asfortranarray([0.1875, 0.0, t_val3])
    end_vals = np.asfortranarray([1.0, s_val3, 1.0])

    nodes = np.asfortranarray([
        [0.1875, 0.0],
        [1.0, 0.0],
        [1.0 - s_val3, s_val3],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (1, 1),
    )
    surface_surface_check(SURFACE1L, SURFACE28Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_29Q():
    s_val1, s_val2 = runtime_utils.real_roots([128, -128, 7])
    t_val1, t_val2 = runtime_utils.real_roots([8, -8, 1])
    start_vals = np.asfortranarray([0.0, 0.0, t_val1, s_val2, 0.0])
    end_vals = np.asfortranarray([1.0, s_val1, t_val2, 1.0, 1.0])

    nodes = np.asfortranarray([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0 - s_val1, s_val1],
        [1.0 - s_val2, s_val2],
        [0.0, 1.0],
    ])
    edge_pairs = (
        (0, 0),
        (0, 1),
        (1, 1),
        (0, 1),
        (0, 2),
    )
    surface_surface_check(SURFACE1L, SURFACE29Q,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces30Q_and_31Q():
    start_vals = np.asfortranarray([13.0 / 56.0, 0.0, 5.0 / 9.0, 0.0])
    end_vals = np.asfortranarray([1.0, 2.0 / 7.0, 1.0, 0.25])

    nodes = np.asfortranarray([
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


def test_surfaces1L_and_7L():
    surface_surface_check_multi(SURFACE1L, SURFACE7L)


def test_surfaces8L_and_29Q():
    surface_surface_check_multi(SURFACE8L, SURFACE29Q)


def test_surfaces1L_and_9L():
    start_vals = np.asfortranarray([0.0, 0.0, 0.0])
    end_vals = np.asfortranarray([1.0, 1.0, 1.0])

    nodes = np.asfortranarray([
        [0.0, 0.125],
        [0.875, 0.0],
        [0.25, 0.75],
    ])
    edge_pairs = (
        (1, 0),
        (1, 1),
        (1, 2),
    )
    surface_surface_check(SURFACE1L, SURFACE9L,
                          start_vals, end_vals, nodes, edge_pairs)


def test_surfaces1L_and_10L():
    surface_surface_check_multi(SURFACE1L, SURFACE10L)


@pytest.mark.parametrize(
    'intersection_info',
    INTERSECTIONS,
    ids=operator.attrgetter('test_id'),
)
def test_intersect(intersection_info):
    id_ = intersection_info.id_
    if intersection_info.note == 'data-unfinished':
        pytest.skip('Intersection does not have all data yet.')

    if id_ in FAILED_CASES_TANGENT:
        kwargs = FAILED_CASES_TANGENT[id_]
        context = check_tangent_manager(**kwargs)
    elif id_ in WIGGLES:
        context = CONFIG.wiggle(WIGGLES[id_])
    else:
        context = runtime_utils.no_op_manager()

    intersected = []
    for curved_poly_info in intersection_info.intersection_infos:
        intersected.append(Intersected(
            curved_poly_info.start_params,
            curved_poly_info.end_params,
            curved_poly_info.intersections,
            curved_poly_info.edge_pairs,
        ))

    surface1 = intersection_info.surface1_info.surface
    surface2 = intersection_info.surface2_info.surface
    with context:
        surface_surface_check_multi(surface1, surface2, *intersected)


if __name__ == '__main__':
    CONFIG.run(globals())
