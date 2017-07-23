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
import contextlib
import itertools

import pytest
import six

import bezier
from bezier import _implicitization
from bezier import curve

import runtime_utils


ALGEBRAIC = curve.IntersectionStrategy.algebraic
GEOMETRIC = curve.IntersectionStrategy.geometric
_, INTERSECTIONS = runtime_utils.surface_intersections_info()
PARALLEL_FAILURE = ('Line segments parallel.',)
BAD_TANGENT = (
    'Curves moving in opposite direction but define '
    'overlapping arcs.')
BAD_TANGENT = (BAD_TANGENT,)
TANGENT_FAILURE = 'The number of candidate intersections is too high.'
WIGGLES = {
    GEOMETRIC: {
        1: 46,  # Established on Ubuntu 16.04
        3: 10,  # Established on AppVeyor (64-bit Python 2.7)
        13: 19,  # Established on Ubuntu 16.04
        22: 37,  # Established on Ubuntu 16.04
        32: 1013,  # Established on Ubuntu 16.04
        33: 1013,  # Established on Ubuntu 16.04
    },
    ALGEBRAIC: {
        3: 12,  # Established on Ubuntu 16.04
        22: 14,  # Established on Ubuntu 16.04
        32: 18,  # Established on Ubuntu 16.04
        33: 18,  # Established on Ubuntu 16.04
    },
}
FAILED_CASES_TANGENT = {
    GEOMETRIC: {
        7: {},
        10: {'parallel': True},
        11: {},
        12: {},
        21: {'bad_tangent': True},
    },
    ALGEBRAIC: {
        6: {},
        7: {},
        10: {},
        11: {},
        12: {},
        15: {},
        17: {},
        21: {},
    },
}
FAILED_CASES_COINCIDENT = {
    GEOMETRIC: {
        4: {},
        5: {'parallel': True},
    },
    ALGEBRAIC: {
        4: {},
        5: {},
    },
}
CONFIG = runtime_utils.Config()


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
def check_tangent_manager(strategy, parallel=False, bad_tangent=False):
    caught_exc = None
    try:
        yield
    except NotImplementedError as exc:
        caught_exc = exc

    assert caught_exc is not None
    exc_args = caught_exc.args
    if strategy is GEOMETRIC:
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
def check_coincident_manager(strategy, parallel=False):
    caught_exc = None
    try:
        yield
    except NotImplementedError as exc:
        caught_exc = exc

    assert caught_exc is not None
    exc_args = caught_exc.args
    if strategy is GEOMETRIC:
        if parallel:
            assert exc_args == PARALLEL_FAILURE
        else:
            assert len(exc_args) == 1
            assert exc_args[0].startswith(TANGENT_FAILURE)
    else:
        assert exc_args == (_implicitization._COINCIDENT_ERR,)


def surface_surface_check(strategy, surface1, surface2, *all_intersected):
    # pylint: disable=too-many-locals
    assert surface1.is_valid
    assert surface2.is_valid

    intersections = surface1.intersect(surface2, strategy=strategy)
    assert len(intersections) == len(all_intersected)
    edges = (
        surface1._get_edges(),
        surface2._get_edges(),
    )

    for intersection, intersected in six.moves.zip(
            intersections, all_intersected):
        assert isinstance(intersection, bezier.CurvedPolygon)

        start_vals = intersected.start_params
        end_vals = intersected.end_params
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
    # pylint: enable=too-many-locals


def _id_func(value):
    if isinstance(value, curve.IntersectionStrategy):
        return 'strategy: {}'.format(value.name)
    else:
        return value.test_id


@pytest.mark.parametrize(
    'strategy,intersection_info',
    itertools.product(
        (GEOMETRIC, ALGEBRAIC),
        INTERSECTIONS,
    ),
    ids=_id_func,
)
def test_intersect(strategy, intersection_info):
    id_ = intersection_info.id_
    if id_ in FAILED_CASES_TANGENT[strategy]:
        kwargs = FAILED_CASES_TANGENT[strategy][id_]
        context = check_tangent_manager(strategy, **kwargs)
    elif id_ in FAILED_CASES_COINCIDENT[strategy]:
        kwargs = FAILED_CASES_COINCIDENT[strategy][id_]
        context = check_coincident_manager(strategy, **kwargs)
    elif id_ in WIGGLES[strategy]:
        context = CONFIG.wiggle(WIGGLES[strategy][id_])
    else:
        context = runtime_utils.no_op_manager()

    surface1 = intersection_info.surface1_info.surface
    surface2 = intersection_info.surface2_info.surface
    with context:
        surface_surface_check(
            strategy, surface1, surface2,
            *intersection_info.intersections)
