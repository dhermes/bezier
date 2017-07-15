# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import

import operator

import pytest
import six

from bezier import _intersection_helpers

import runtime_utils


CONFIG = runtime_utils.Config()
S_PROP = operator.attrgetter('s')
CURVES, INTERSECTIONS = runtime_utils.curve_intersections_info()
WIGGLES = {
    ('8', '27'): 41,  # Established on Ubuntu 16.04
    ('14', '16'): 4,  # Established on Ubuntu 16.04 (Less than 8)
    ('20', '21'): 7,  # Established on Ubuntu 16.04 (Less than 8)
    ('21', '22'): 11,  # Established on Ubuntu 16.04
    ('32', '33'): 3,  # Established on Ubuntu 16.04 (Less than 8)
    ('50', '54'): 91,  # Established on Ubuntu 16.04
    ('51', '54'): 1013,  # Established on Ubuntu 16.04
    ('52', '54'): 91,  # Established on Ubuntu 16.04
    ('53', '54'): 1013,  # Established on Ubuntu 16.04
}
FAILURE_NOT_IMPLEMENTED = (
    ('1', '24'),  # The number of candidate intersections is too high. (24)
    ('14', '15'),  # Line segments parallel.
    ('28', '29'),  # The number of candidate intersections is too high. (22)
    ('56', '57'),  # The number of candidate intersections is too high. (20)
)
INCORRECT_COUNT = (
    ('38', '39'),
)


class IncorrectCount(ValueError):
    """Custom exception for a "very bad" answer.

    This should be raised when the **computed** number of intersections
    disagrees with the actual number of intersections.
    """


def _get_params(intersection_info):
    s_vals = intersection_info['curve1_params']
    num_s, = s_vals.shape

    t_vals = intersection_info['curve2_params']
    assert t_vals.shape == (num_s,)

    intersection_pts = intersection_info['intersections']
    if num_s == 0:
        assert intersection_pts.size == 0
    else:
        assert intersection_pts.shape == (num_s, 2)

    return s_vals, t_vals, intersection_pts


def _intersection_curves(num_s, curve_id1, curve_id2, curve1, curve2):
    strategy = _intersection_helpers.IntersectionStrategy.geometric
    intersections = _intersection_helpers.all_intersections(
        [(curve1, curve2)], strategy=strategy)

    # Make we have the right number of intersections.
    if len(intersections) != num_s:
        raise IncorrectCount(
            'Received wrong number of intersections',
            len(intersections), 'Expected', num_s,
            'Curve 1', curve_id1, 'Curve 2', curve_id2)

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


def _intersections_check(curve_id1, curve_id2, intersection_info):
    # Get curve instances.
    curve1 = CURVES[curve_id1].curve
    curve2 = CURVES[curve_id2].curve

    # Get / verify info for the intersections.
    s_vals, t_vals, intersection_pts = _get_params(intersection_info)

    # Actually intersect the curves.
    intersections = _intersection_curves(
        s_vals.size, curve_id1, curve_id2, curve1, curve2)

    # Check that each intersection is as expected.
    info = six.moves.zip(intersections, s_vals, t_vals, intersection_pts)
    for info_tuple in info:
        _intersection_check(info_tuple, curve1, curve2)


@pytest.mark.parametrize(
    'intersection_info',
    INTERSECTIONS,
    ids=runtime_utils.curve_id_func,
)
def test_intersect(intersection_info):
    curve_id1 = intersection_info['curve1']
    curve_id2 = intersection_info['curve2']
    id_pair = (curve_id1, curve_id2)

    if id_pair in FAILURE_NOT_IMPLEMENTED:
        assert intersection_info['type'] in ('tangent', 'coincident')
        context = pytest.raises(NotImplementedError)
    elif id_pair in INCORRECT_COUNT:
        assert intersection_info['type'] == 'tangent'
        context = pytest.raises(IncorrectCount)
    elif id_pair in WIGGLES:
        context = CONFIG.wiggle(WIGGLES[id_pair])
    else:
        context = runtime_utils.no_op_manager()

    with context:
        _intersections_check(curve_id1, curve_id2, intersection_info)
