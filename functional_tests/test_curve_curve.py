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

import numpy as np
import pytest
import six

import bezier.curve
from bezier import _intersection_helpers

import runtime_utils
from runtime_utils import CurveIntersectionType


SPACING = np.spacing  # pylint: disable=no-member
CONFIG = runtime_utils.Config()
GEOMETRIC = bezier.curve.IntersectionStrategy.geometric
S_PROP = operator.attrgetter('s')
_, INTERSECTIONS = runtime_utils.curve_intersections_info()
ULPS_ALLOWED = 3.0
# NOTE: Set a threshold for values that are "approximately" zero.
#       This is an **absolute** error rather than a relative error
#       since relative to zero error is always infinite.
ZERO_THRESHOLD = 0.5**52
ZERO_MISS = object()
ULPS_ALLOWED_OVERRIDE = {
    GEOMETRIC: {
        12: {
            (0, 1): 4,  # Established on Ubuntu 16.04
        },
        17: {
            (0, 7): ZERO_MISS,  # Established on Ubuntu 16.04
            (2, 1): 7,  # Established on Ubuntu 16.04
            (3, 1): 4,  # Established on Ubuntu 16.04
            (3, 2): 7,  # Established on Ubuntu 16.04
        },
        18: {
            (0, 6): 5,  # Established on Ubuntu 16.04
            (1, 0): 8,  # Established on Ubuntu 16.04
            (1, 2): 10,  # Established on Ubuntu 16.04
            (1, 3): 11,  # Established on Ubuntu 16.04
            (1, 6): 4,  # Established on Ubuntu 16.04
            (2, 2): 4,  # Established on Ubuntu 16.04
            (3, 3): ZERO_MISS,  # Established on Ubuntu 16.04
            (3, 5): ZERO_MISS,  # Established on Ubuntu 16.04
            (3, 7): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        23: {
            (0, 0): 14,  # Established on Ubuntu 16.04
            (0, 1): 41,  # Established on Ubuntu 16.04
            (0, 2): 14,  # Established on Ubuntu 16.04
            (1, 0): 16,  # Established on Ubuntu 16.04
            (1, 1): 21,  # Established on Ubuntu 16.04
            (1, 2): 16,  # Established on Ubuntu 16.04
        },
        25: {
            (0, 6): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        37: {
            (0, 1): 91,  # Established on Ubuntu 16.04
        },
        38: {
            (0, 1): 1013,  # Established on Ubuntu 16.04
        },
        39: {
            (0, 1): 91,  # Established on Ubuntu 16.04
        },
        40: {
            (0, 1): 1013,  # Established on Ubuntu 16.04
        },
    },
}
TANGENT_SUCCESS = {
    GEOMETRIC: (
        4,
        14,
        19,
        41,
    ),
}
COINCIDENT_SUCCESS = {
    GEOMETRIC: (
        33,
        34,
        35,
    ),
}
INCORRECT_COUNT = {
    GEOMETRIC: (
        31,
    ),
}


class IncorrectCount(ValueError):
    """Custom exception for a "very bad" answer.

    This should be raised when the **computed** number of intersections
    disagrees with the actual number of intersections.
    """


def get_sorted_intersections(intersection_info, strategy):
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


def intersection_values(intersection_info, strategy):
    # Actually intersect the curves.
    intersections = get_sorted_intersections(intersection_info, strategy)

    curve1 = intersection_info.curve1
    curve2 = intersection_info.curve2

    # Put the computed and expected values into matrices to be compared.
    s_vals, t_vals, intersection_pts = intersection_info.params
    values_shape = (intersection_info.num_params, 8)
    computed = np.zeros(values_shape, order='F')
    exact = np.zeros(values_shape, order='F')

    info = six.moves.zip(intersections, s_vals, t_vals, intersection_pts)
    for index, (intersection, s_val, t_val, point) in enumerate(info):
        assert intersection.first is curve1
        assert intersection.second is curve2

        computed[index, (0, 1)] = intersection.s, intersection.t
        exact[index, (0, 1)] = s_val, t_val

        computed_point = intersection.get_point()
        computed[index, (2, 3)] = computed_point[0, :]
        exact[index, (2, 3)] = point

        point_on1 = curve1.evaluate(s_val)
        computed[index, (4, 5)] = point_on1[0, :]
        exact[index, (4, 5)] = point

        point_on2 = curve2.evaluate(t_val)
        computed[index, (6, 7)] = point_on2[0, :]
        exact[index, (6, 7)] = point

    return computed, exact


def error_multipliers(intersection_info, shape, strategy):
    zero_misses = []
    multipliers = ULPS_ALLOWED * np.ones(shape, order='F')
    override = ULPS_ALLOWED_OVERRIDE[strategy].get(intersection_info.id_)
    if override is not None:
        for index_tuple, value in six.iteritems(override):
            if value is ZERO_MISS:
                multipliers[index_tuple] = np.inf
                zero_misses.append(index_tuple)
            else:
                multipliers[index_tuple] = value

    return zero_misses, multipliers


def check_intersect(intersection_info, strategy):
    computed, exact = intersection_values(intersection_info, strategy)
    zero_misses, multipliers = error_multipliers(
        intersection_info, exact.shape, strategy)

    # Make sure our errors fall under the number of "allowed" ULPs.
    allowed_errors = np.abs(multipliers * SPACING(exact))
    errors = np.abs(exact - computed)
    assert np.all(errors <= allowed_errors)

    # If there are any ``exact`` zeros that have been missed, check
    # that we fall under the **absolute** threshold for them.
    if zero_misses:
        for index_tuple in zero_misses:
            assert exact[index_tuple] == 0.0
            assert np.abs(computed[index_tuple]) < ZERO_THRESHOLD


def check_no_intersect(intersection_info, strategy):
    computed, exact = intersection_values(intersection_info, strategy)
    assert computed.size == 0
    assert exact.size == 0


def check_tangent(intersection_info, strategy):
    id_ = intersection_info.id_
    if id_ in TANGENT_SUCCESS[strategy]:
        check_intersect(intersection_info, strategy)
    else:
        with pytest.raises(NotImplementedError):
            intersection_values(intersection_info, strategy)


def check_coincident(intersection_info, strategy):
    id_ = intersection_info.id_
    if id_ in COINCIDENT_SUCCESS[strategy]:
        check_intersect(intersection_info, strategy)
    else:
        with pytest.raises(NotImplementedError):
            intersection_values(intersection_info, strategy)


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
    # Actually try to intersect the curves.
    intersection_type = intersection_info.type_

    if id_ in INCORRECT_COUNT[strategy]:
        assert intersection_info.type_ == CurveIntersectionType.tangent
        with pytest.raises(IncorrectCount):
            check_intersect(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.tangent:
        check_tangent(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.coincident:
        check_coincident(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.standard:
        check_intersect(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.no_intersection:
        check_no_intersect(intersection_info, strategy)
    else:
        raise ValueError(
            'Unexpected intersection type', intersection_type)
