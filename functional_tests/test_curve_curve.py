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

import itertools
import operator

import pytest
import six

import bezier.curve
from bezier import _intersection_helpers

import runtime_utils
from runtime_utils import CurveIntersectionType


CONFIG = runtime_utils.Config()
GEOMETRIC = bezier.curve.IntersectionStrategy.geometric
S_PROP = operator.attrgetter('s')
_, INTERSECTIONS = runtime_utils.curve_intersections_info()
WIGGLES = {
    12: 4,  # Established on Ubuntu 16.04 (Less than 8)
    17: 7,  # Established on Ubuntu 16.04 (Less than 8)
    18: 11,  # Established on Ubuntu 16.04
    23: 41,  # Established on Ubuntu 16.04
    28: 3,  # Established on Ubuntu 16.04 (Less than 8)
    37: 91,  # Established on Ubuntu 16.04
    38: 1013,  # Established on Ubuntu 16.04
    39: 91,  # Established on Ubuntu 16.04
    40: 1013,  # Established on Ubuntu 16.04
}
FAILURE_NOT_IMPLEMENTED = (
    11,  # Line segments parallel.
    20,  # The number of candidate intersections is too high. (24)
    24,  # The number of candidate intersections is too high. (22)
    42,  # The number of candidate intersections is too high. (20)
)
NOT_IMPLEMENTED_TYPES = (
    CurveIntersectionType.tangent,
    CurveIntersectionType.coincident,
)
INCORRECT_COUNT = (
    31,
)


class IncorrectCount(ValueError):
    """Custom exception for a "very bad" answer.

    This should be raised when the **computed** number of intersections
    disagrees with the actual number of intersections.
    """


def _intersection_curves(intersection_info):
    strategy = _intersection_helpers.IntersectionStrategy.geometric
    curve1 = intersection_info.curve1
    curve2 = intersection_info.curve2
    intersections = _intersection_helpers.all_intersections(
        [(curve1, curve2)], strategy=strategy)

    # Make we have the right number of intersections.
    if len(intersections) != intersection_info.num_params:
        raise IncorrectCount(
            'Received wrong number of intersections',
            len(intersections), 'Expected', intersection_info.num_params,
            intersection_info.test_id)

    # Sort the intersections by s-value.
    intersections.sort(key=S_PROP)
    return intersections


def _intersection_check(info_tuple, curve1, curve2):
    intersection, s_val, t_val, point = info_tuple
    assert intersection.first is curve1
    assert intersection.second is curve2

    CONFIG.assert_close(intersection.s, s_val)
    CONFIG.assert_close(intersection.t, t_val)

    computed_point = intersection.get_point()
    CONFIG.assert_close(computed_point[0, 0], point[0])
    CONFIG.assert_close(computed_point[0, 1], point[1])

    point_on1 = curve1.evaluate(s_val)
    CONFIG.assert_close(point_on1[0, 0], point[0])
    CONFIG.assert_close(point_on1[0, 1], point[1])

    point_on2 = curve2.evaluate(t_val)
    CONFIG.assert_close(point_on2[0, 0], point[0])
    CONFIG.assert_close(point_on2[0, 1], point[1])


def _intersections_check(intersection_info):
    # Actually intersect the curves.
    intersections = _intersection_curves(intersection_info)

    # Check that each intersection is as expected.
    s_vals, t_vals, intersection_pts = intersection_info.params
    info = six.moves.zip(intersections, s_vals, t_vals, intersection_pts)
    curve1 = intersection_info.curve1
    curve2 = intersection_info.curve2
    for info_tuple in info:
        _intersection_check(info_tuple, curve1, curve2)


@pytest.mark.parametrize(
    'strategy,intersection_info',
    itertools.product(
        (GEOMETRIC,),
        INTERSECTIONS,
    ),
    ids=runtime_utils.id_func,
)
def test_intersect(strategy, intersection_info):
    id_ = intersection_info.id_
    if id_ in FAILURE_NOT_IMPLEMENTED:
        assert intersection_info.type_ in NOT_IMPLEMENTED_TYPES
        context = pytest.raises(NotImplementedError)
    elif id_ in INCORRECT_COUNT:
        assert intersection_info.type_ == CurveIntersectionType.tangent
        context = pytest.raises(IncorrectCount)
    elif id_ in WIGGLES:
        context = CONFIG.wiggle(WIGGLES[id_])
    else:
        context = runtime_utils.no_op_manager()

    with context:
        _intersections_check(intersection_info)
