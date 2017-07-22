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
    7: {},
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


def curved_polygon_edges(intersection, edges):
    edges1, edges2 = edges
    all_edges = edges1 + edges2
    # Re-sort the edges to be in the same order independent of strategy.
    edge_list = intersection._edges
    edge_info = [(all_edges.index(edge.root), edge.start, edge.end)
                 for edge in edge_list]
    index = edge_info.index(min(edge_info))
    return edge_list[index:] + edge_list[:index]


@contextlib.contextmanager
def check_tangent_manager(parallel=False, bad_tangent=False):
    caught_exc = None
    try:
        yield
    except NotImplementedError as exc:
        caught_exc = exc

    assert caught_exc is not None
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
def check_coincident_manager(parallel=False):
    caught_exc = None
    try:
        yield
    except NotImplementedError as exc:
        caught_exc = exc

    assert caught_exc is not None
    exc_args = caught_exc.args
    if STRATEGY is GEOMETRIC:
        if parallel:
            assert exc_args == PARALLEL_FAILURE
        else:
            assert len(exc_args) == 1
            assert exc_args[0].startswith(TANGENT_FAILURE)
    else:
        assert exc_args == (_implicitization._COINCIDENT_ERR,)


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


@pytest.mark.parametrize(
    'intersection_info',
    INTERSECTIONS,
    ids=operator.attrgetter('test_id'),
)
def test_intersect(intersection_info):
    id_ = intersection_info.id_
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

    intersected = [
        Intersected(
            curved_poly_info.start_params,
            curved_poly_info.end_params,
            curved_poly_info.intersections,
            curved_poly_info.edge_pairs,
        )
        for curved_poly_info in intersection_info.intersection_infos
    ]

    surface1 = intersection_info.surface1_info.surface
    surface2 = intersection_info.surface2_info.surface
    with context:
        surface_surface_check_multi(surface1, surface2, *intersected)


if __name__ == '__main__':
    CONFIG.run(globals())
