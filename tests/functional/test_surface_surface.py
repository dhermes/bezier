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

import numpy as np
import pytest
import six

import bezier
from bezier import _algebraic_intersection
from bezier import _geometric_intersection
from bezier import _surface_helpers
from bezier import curve
from tests.functional import utils


ALGEBRAIC = curve.IntersectionStrategy.ALGEBRAIC
GEOMETRIC = curve.IntersectionStrategy.GEOMETRIC
_, INTERSECTIONS = utils.surface_intersections_info()
PARALLEL_FAILURE = (_geometric_intersection._SEGMENTS_PARALLEL,)
BAD_TANGENT = (_surface_helpers._BAD_TANGENT,)
TOO_MANY = _geometric_intersection._TOO_MANY_TEMPLATE
WIGGLES = {
    GEOMETRIC: {
        1: 48,  # Established on CentOS 5 (i686 Docker image)
        3: 10,  # Established on AppVeyor (64-bit Python 2.7)
        8: 18,  # Established on CentOS 5 (i686 Docker image)
        13: 19,  # Established on Ubuntu 16.04
        22: 37,  # Established on Ubuntu 16.04
        32: 1013,  # Established on Ubuntu 16.04
        33: 1013,  # Established on Ubuntu 16.04
        34: 19,  # Established on Ubuntu 16.04 on CircleCI
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
        7: {'too_many': 66},
        10: {'parallel': True},
        11: {'too_many': 78},
        12: {'too_many': 76},
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
        4: {'too_many': 96},
        5: {'parallel': True},
    },
    ALGEBRAIC: {
        4: {},
        5: {},
    },
}
CONFIG = utils.Config()


def find_edge_index(edge, root_edges):
    for index, root_edge in enumerate(root_edges):
        if root_edge.degree != edge.degree:
            continue
        specialized = root_edge.specialize(edge.start, edge.end)
        if np.allclose(specialized._nodes, edge._nodes):
            return index

    raise RuntimeError('No match found.')


def curved_polygon_edges(intersection, root_edges_pair):
    edges1, edges2 = root_edges_pair
    root_edges = edges1 + edges2
    # Re-sort the edges to be in the same order independent of strategy.
    edge_list = intersection._edges
    edge_info = [
        (find_edge_index(edge, root_edges), edge.start, edge.end)
        for edge in edge_list
    ]
    index = edge_info.index(min(edge_info))
    return edge_list[index:] + edge_list[:index]


@contextlib.contextmanager
def check_tangent_manager(
        strategy, parallel=False, bad_tangent=False, too_many=None):
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
            err_msg = TOO_MANY.format(too_many)
            assert exc_args == (err_msg,)
    else:
        assert len(exc_args) == 2
        assert exc_args[0] == _algebraic_intersection._NON_SIMPLE_ERR


@contextlib.contextmanager
def check_coincident_manager(strategy, parallel=False, too_many=None):
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
            err_msg = TOO_MANY.format(too_many)
            assert exc_args == (err_msg,)
    else:
        assert exc_args == (_algebraic_intersection._COINCIDENT_ERR,)


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
        if isinstance(intersection, bezier.CurvedPolygon):
            if intersected.parent.id_ in (18, 19, 30):
                # NOTE: This branch is a **temporary** hack. Once resolved,
                #       the properties (such as ``edge_nodes``) can be removed
                #       from ``SurfaceIntersectionInfo``.
                assert isinstance(intersected, utils.SurfaceIntersectionInfo)
            else:
                assert isinstance(intersected, utils.CurvedPolygonInfo)
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
                expected_root = edges[surf_index][edge_index]
                specialized = expected_root.specialize(start_val, end_val)
                assert np.allclose(specialized._nodes, edge._nodes)

                CONFIG.assert_close(edge.start, start_val)
                CONFIG.assert_close(edge.end, end_val)
                CONFIG.assert_close(edge._nodes[0, 0], node[0])
                CONFIG.assert_close(edge._nodes[0, 1], node[1])
        else:
            assert isinstance(intersection, bezier.Surface)
            assert isinstance(intersected, utils.SurfaceIntersectionInfo)
            if intersected.first:
                assert intersection is surface1
            else:
                assert intersection is surface2
    # pylint: enable=too-many-locals


@pytest.mark.parametrize(
    'strategy,intersection_info',
    itertools.product(
        (GEOMETRIC, ALGEBRAIC),
        INTERSECTIONS,
    ),
    ids=utils.id_func,
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
        context = utils.no_op_manager()

    surface1 = intersection_info.surface1
    surface2 = intersection_info.surface2
    with context:
        surface_surface_check(
            strategy, surface1, surface2,
            *intersection_info.intersections)
