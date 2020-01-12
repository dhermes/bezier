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

import contextlib
import itertools

import numpy as np
import pytest

import bezier
from bezier import _algebraic_intersection
from bezier import _py_geometric_intersection
from bezier import _py_triangle_helpers
from bezier import _py_triangle_intersection
from bezier import curve
from tests import utils as base_utils
from tests.functional import utils

ALGEBRAIC = curve.IntersectionStrategy.ALGEBRAIC
GEOMETRIC = curve.IntersectionStrategy.GEOMETRIC
_, INTERSECTIONS = utils.triangle_intersections_info()
SAME_CURVATURE = (_py_triangle_helpers._SAME_CURVATURE,)
TOO_MANY = _py_geometric_intersection._TOO_MANY_TEMPLATE
WIGGLES = {
    GEOMETRIC: {
        1: 48,  # Established on CentOS 5 (i686 Docker image)
        3: 10,  # Established on AppVeyor (64-bit Python 2.7)
        8: 18,  # Established on CentOS 5 (i686 Docker image)
        13: 19,  # Established on Ubuntu 16.04
        21: 15,  # Established on Ubuntu 16.04
        22: 37,  # Established on Ubuntu 16.04
        32: 1013,  # Established on Ubuntu 16.04
        33: 1013,  # Established on Ubuntu 16.04
        34: 19,  # Established on Ubuntu 16.04 on CircleCI
        49: 86670,  # Established on Ubuntu 16.04
        50: 15222,  # Established on Ubuntu 16.04
        51: 285086,  # Established on CentOS 5 (i686 Docker image)
    },
    ALGEBRAIC: {
        3: 12,  # Established on Ubuntu 16.04
        8: 20,  # Established on Ubuntu 18.04
        22: 14,  # Established on Ubuntu 16.04
        32: 18,  # Established on Ubuntu 16.04
        33: 18,  # Established on Ubuntu 16.04
        36: 9,  # Established on Ubuntu 16.04 on CircleCI and macOS.
        49: 120841,  # Established on CentOS 5 (i686 Docker image)
        50: 86670,  # Established on Ubuntu 16.04
    },
}
FAILED_CASES_TANGENT = {
    GEOMETRIC: {},
    ALGEBRAIC: {
        6: {},
        7: {},
        10: {},
        11: {},
        12: {},
        15: {},
        17: {},
        21: {},
        37: {},
        41: {},
        42: {},
        48: {},
    },
}
FAILED_CASES_COINCIDENT = {
    GEOMETRIC: {},
    ALGEBRAIC: {4: {}, 5: {}, 43: {}, 44: {}, 45: {}, 46: {}, 47: {}, 51: {}},
}
FAILED_CASES_BAD_EDGES = {GEOMETRIC: (), ALGEBRAIC: ()}
FAILED_CASES_PRECISION = {GEOMETRIC: (), ALGEBRAIC: ()}
FAILED_CASES_BAD_DUPLICATE = {
    GEOMETRIC: (),
    ALGEBRAIC: (53, 54, 55, 56, 57, 59, 60, 62, 63, 64, 66, 67, 69, 70),
}
FAILED_CASES_CONSECUTIVE_SEGMENTS = {
    GEOMETRIC: (55, 56, 59, 61, 63, 65, 66, 67, 70, 71, 72),
    ALGEBRAIC: (58, 61, 65, 71),
}
INCORRECT_COUNT = {GEOMETRIC: (52, 62, 64, 69), ALGEBRAIC: (68,)}
if base_utils.IS_LINUX and not base_utils.IS_64_BIT:
    FAILED_CASES_PRECISION[GEOMETRIC] += (54, 57)
    FAILED_CASES_CONSECUTIVE_SEGMENTS[GEOMETRIC] += (68,)
    INCORRECT_COUNT[GEOMETRIC] += (53, 58, 60)
    INCORRECT_COUNT[ALGEBRAIC] += (52, 72)
else:
    FAILED_CASES_PRECISION[GEOMETRIC] += (68,)
    FAILED_CASES_BAD_EDGES[ALGEBRAIC] += (52,)
    FAILED_CASES_CONSECUTIVE_SEGMENTS[GEOMETRIC] += (53, 54, 57, 58, 60)
    FAILED_CASES_CONSECUTIVE_SEGMENTS[ALGEBRAIC] += (72,)

CONFIG = utils.Config()


def curved_polygon_edges(intersection):
    # Re-sort the edges to be in the same order independent of strategy.
    edge_info = intersection._metadata
    index = edge_info.index(min(edge_info))
    edge_info = edge_info[index:] + edge_info[:index]
    edge_list = intersection._edges
    edge_list = edge_list[index:] + edge_list[:index]
    return edge_info, edge_list


@contextlib.contextmanager
def check_tangent_manager(strategy, too_many=None):
    caught_exc = None
    try:
        yield

    except NotImplementedError as exc:
        caught_exc = exc
    assert caught_exc is not None
    exc_args = caught_exc.args
    if strategy is GEOMETRIC:
        err_msg = TOO_MANY.format(too_many)
        assert exc_args == (err_msg,)
    else:
        assert len(exc_args) == 2
        assert exc_args[0] == _algebraic_intersection._NON_SIMPLE_ERR


@contextlib.contextmanager
def check_coincident_manager(strategy, too_many=None, curvature=False):
    caught_exc = None
    try:
        yield

    except NotImplementedError as exc:
        caught_exc = exc
    assert caught_exc is not None
    exc_args = caught_exc.args
    if strategy is GEOMETRIC:
        if curvature:
            assert exc_args == SAME_CURVATURE
        else:
            err_msg = TOO_MANY.format(too_many)
            assert exc_args == (err_msg,)
    else:
        assert exc_args == (_algebraic_intersection._COINCIDENT_ERR,)


@contextlib.contextmanager
def check_bad_edges_manager():
    caught_exc = None
    try:
        yield

    except RuntimeError as exc:
        caught_exc = exc
    assert caught_exc is not None
    exc_args = caught_exc.args
    assert exc_args[0] == "Unexpected number of edges"


@contextlib.contextmanager
def check_precision_manager():
    caught_exc = None
    try:
        yield

    except AssertionError as exc:
        caught_exc = exc

    assert caught_exc is not None


@contextlib.contextmanager
def check_duplicate_manager():
    caught_exc = None
    try:
        yield

    except ValueError as exc:
        caught_exc = exc

    assert caught_exc is not None
    exc_args = caught_exc.args
    assert exc_args[0] == "Duplicate not among uniques"


@contextlib.contextmanager
def check_consecutive_manager():
    caught_exc = None
    try:
        yield

    except ValueError as exc:
        caught_exc = exc

    assert caught_exc is not None
    exc_args = caught_exc.args
    assert exc_args[0] == _py_triangle_intersection.SEGMENTS_SAME_EDGE
    assert len(exc_args) == 3


def extra_verify(strategy, intersections):
    """Do extra verification on a list of intersections.

    This is intended to be used for cases when the "regular" Python
    verification was not run, e.g. if the Fortran speedups were used.

    Args:
        strategy (.IntersectionStrategy): The strategy that was used to
            intersect edges.
        intersections (List[Union[~bezier.curved_polygon.CurvedPolygon, \
            ~bezier.triangle.Triangle]]): List of intersections (possibly empty).
    """
    if strategy == GEOMETRIC and bezier._HAS_SPEEDUP:
        edge_infos = [
            curved_polygon._metadata
            for curved_polygon in intersections
            if isinstance(curved_polygon, bezier.CurvedPolygon)
        ]
        _py_triangle_intersection.verify_edge_segments(edge_infos)


def triangle_triangle_check(strategy, triangle1, triangle2, *all_intersected):
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-locals
    assert triangle1.is_valid
    assert triangle2.is_valid
    intersections = triangle1.intersect(triangle2, strategy=strategy)
    extra_verify(strategy, intersections)

    if len(intersections) != len(all_intersected):
        raise utils.IncorrectCount(
            "Received wrong number of intersections",
            len(intersections),
            "Expected",
            len(all_intersected),
        )

    all_edges = triangle1._get_edges() + triangle2._get_edges()
    for intersection, intersected in zip(intersections, all_intersected):
        if isinstance(intersection, bezier.CurvedPolygon):
            assert isinstance(intersected, utils.CurvedPolygonInfo)
            edge_info, int_edges = curved_polygon_edges(intersection)
            num_edges = len(edge_info)
            if num_edges != len(intersected.edge_list):
                raise utils.IncorrectCount(
                    "Received wrong number of edges within an intersection",
                    num_edges,
                    "Expected",
                    len(intersected.edge_list),
                )

            assert num_edges == len(intersected.start_params)
            assert num_edges == len(intersected.end_params)
            assert num_edges == len(intersected.nodes)
            for index in range(num_edges):
                edge_triple = edge_info[index]
                edge_index = intersected.edge_list[index]
                assert edge_triple[0] == edge_index
                start_val = intersected.start_params[index]
                end_val = intersected.end_params[index]
                CONFIG.assert_close(edge_triple[1], start_val)
                CONFIG.assert_close(edge_triple[2], end_val)
                expected_root = all_edges[edge_index]
                specialized = expected_root.specialize(start_val, end_val)
                edge = int_edges[index]
                assert np.allclose(specialized._nodes, edge._nodes)
                node = intersected.nodes[index, :]
                CONFIG.assert_close(edge._nodes[0, 0], node[0])
                CONFIG.assert_close(edge._nodes[1, 0], node[1])
        else:
            assert isinstance(intersection, bezier.Triangle)
            assert isinstance(intersected, utils.TriangleIntersectionInfo)
            if intersected.first:
                assert intersection is triangle1
            else:
                assert intersection is triangle2


@pytest.mark.parametrize(
    "strategy,intersection_info",
    itertools.product((GEOMETRIC, ALGEBRAIC), INTERSECTIONS),
    ids=utils.id_func,
)
def test_intersect(strategy, intersection_info):
    id_ = intersection_info.id_
    if id_ in INCORRECT_COUNT[strategy]:
        context = pytest.raises(utils.IncorrectCount)
    elif id_ in FAILED_CASES_TANGENT[strategy]:
        kwargs = FAILED_CASES_TANGENT[strategy][id_]
        context = check_tangent_manager(strategy, **kwargs)
    elif id_ in FAILED_CASES_COINCIDENT[strategy]:
        kwargs = FAILED_CASES_COINCIDENT[strategy][id_]
        context = check_coincident_manager(strategy, **kwargs)
    elif id_ in FAILED_CASES_BAD_EDGES[strategy]:
        context = check_bad_edges_manager()
    elif id_ in FAILED_CASES_PRECISION[strategy]:
        context = check_precision_manager()
    elif id_ in FAILED_CASES_BAD_DUPLICATE[strategy]:
        context = check_duplicate_manager()
    elif id_ in FAILED_CASES_CONSECUTIVE_SEGMENTS[strategy]:
        context = check_consecutive_manager()
    elif id_ in WIGGLES[strategy]:
        context = CONFIG.wiggle(WIGGLES[strategy][id_])
    else:
        context = utils.no_op_manager()
    triangle1 = intersection_info.triangle1
    triangle2 = intersection_info.triangle2
    with context:
        triangle_triangle_check(
            strategy, triangle1, triangle2, *intersection_info.intersections
        )
