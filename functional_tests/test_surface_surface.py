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
    3: 9,
    13: 19,
    22: 37,
    32: 1013,
    33: 1013,
}
FAILED_CASES_TANGENT = {
    10: {'parallel': True},
    11: {},
    12: {},
    21: {'bad_tangent': True},
}
FAILED_CASES_COINCIDENT = {
    4: {},
    5: {'parallel': True},
}
CONFIG = runtime_utils.Config()

SURFACE1L = SURFACES['1L'].surface
SURFACE3L = SURFACES['3L'].surface
SURFACE5L = SURFACES['5L'].surface
SURFACE6L = SURFACES['6L'].surface
SURFACE9L = SURFACES['9L'].surface
SURFACE1Q = SURFACES['1Q'].surface
SURFACE3Q = SURFACES['3Q'].surface
SURFACE4Q = SURFACES['4Q'].surface
SURFACE5Q = SURFACES['5Q'].surface
SURFACE10Q = SURFACES['10Q'].surface
SURFACE14Q = SURFACES['14Q'].surface
SURFACE17Q = SURFACES['17Q'].surface
SURFACE19Q = SURFACES['19Q'].surface


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


def check_coincident(caught_exc, parallel=False):
    exc_args = caught_exc.args
    if STRATEGY is GEOMETRIC:
        if parallel:
            assert exc_args == PARALLEL_FAILURE
        else:
            assert len(exc_args) == 1
            assert exc_args[0].startswith(TANGENT_FAILURE)
    else:
        assert exc_args == (_implicitization._COINCIDENT_ERR,)


@contextlib.contextmanager
def check_coincident_manager(**kwargs):
    caught_exc = None
    try:
        yield
    except NotImplementedError as exc:
        caught_exc = exc

    assert caught_exc is not None
    check_coincident(caught_exc, **kwargs)


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
    elif id_ in FAILED_CASES_COINCIDENT:
        kwargs = FAILED_CASES_COINCIDENT[id_]
        context = check_coincident_manager(**kwargs)
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
