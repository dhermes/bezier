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

from bezier import _implicitization
from bezier import _intersection_helpers
import bezier.curve

import runtime_utils
from runtime_utils import CurveIntersectionType


SPACING = np.spacing  # pylint: disable=no-member
CONFIG = runtime_utils.Config()
GEOMETRIC = bezier.curve.IntersectionStrategy.geometric
ALGEBRAIC = bezier.curve.IntersectionStrategy.algebraic
S_PROP = operator.attrgetter('s')
_, INTERSECTIONS = runtime_utils.curve_intersections_info()
INCORRECT_VALUES_MSG = """\
Multipliers were:
{}

However, the actual observed (ULP relative) errors were:
{}

The issue seems to be in the following index pair(s):
{}
"""
ULPS_ALLOWED = 3.0
# NOTE: Set a threshold for values that are "approximately" zero.
#       This is an **absolute** error rather than a relative error
#       since relative to zero error is always infinite.
ZERO_THRESHOLD = 0.5**52
ZERO_MISS = object()
# NOTE: We use units of least precision (ULP) as error. These
#       are for the very rare cases where the computed values
#       differ from the actual values by more than 3 ULPs.
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
    ALGEBRAIC: {
        2: {
            (1, 3): 4,  # Established on Ubuntu 16.04
        },
        17: {
            (0, 7): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        18: {
            (0, 6): 5,  # Established on Ubuntu 16.04
            (1, 6): 4,  # Established on Ubuntu 16.04
            (3, 3): ZERO_MISS,  # Established on Ubuntu 16.04
            (3, 5): ZERO_MISS,  # Established on Ubuntu 16.04
            (3, 7): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        22: {
            (0, 0): 12,  # Established on Ubuntu 16.04
            (0, 1): 29,  # Established on Ubuntu 16.04
            (0, 2): 9,  # Established on Ubuntu 16.04
        },
        23: {
            (0, 1): 4,  # Established on Ubuntu 16.04
        },
        25: {
            (0, 6): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        37: {
            (0, 1): 165,  # Established on Ubuntu 16.04
        },
        38: {
            (0, 1): 18,  # Established on Ubuntu 16.04
        },
        39: {
            (0, 1): 165,  # Established on Ubuntu 16.04
        },
        40: {
            (0, 1): 18,  # Established on Ubuntu 16.04
        },
    },
}
NON_SIMPLE_ERR = _implicitization._NON_SIMPLE_ERR
SEGMENTS_PARALLEL = (_intersection_helpers._SEGMENTS_PARALLEL,)
TOO_MANY = _intersection_helpers._TOO_MANY_TEMPLATE
COINCIDENT_ERR = (_implicitization._COINCIDENT_ERR,)
TANGENT_OVERRIDES = {
    GEOMETRIC: {
        4: {'success': True},
        11: {'parallel': True},
        14: {'success': True},
        19: {'success': True},
        24: {'too_many': 22},
        41: {'success': True},
        42: {'too_many': 20},
    },
    ALGEBRAIC: {},
}
COINCIDENT_OVERRIDES = {
    GEOMETRIC: {
        20: {'too_many': 24},
        33: {'success': True},
        34: {'success': True},
        35: {'success': True},
    },
    ALGEBRAIC: {},
}
INCORRECT_COUNT = {
    GEOMETRIC: (
        31,
    ),
    ALGEBRAIC: {},
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

    # Put the computed and expected values into matrices to be compared.
    s_vals, t_vals, intersection_pts = intersection_info.params
    computed = np.zeros((intersection_info.num_params, 8), order='F')
    exact = np.zeros(computed.shape, order='F')

    info = six.moves.zip(intersections, s_vals, t_vals, intersection_pts)
    for index, (intersection, s_val, t_val, point) in enumerate(info):
        assert intersection.first is intersection_info.curve1
        assert intersection.second is intersection_info.curve2

        computed[index, (0, 1)] = intersection.s, intersection.t
        exact[index, (0, 1)] = s_val, t_val

        # Make sure the "best" point is close to the exact one.
        computed[index, (2, 3)] = intersection.get_point()
        exact[index, (2, 3)] = point

        # Make sure the point corresponding to the parameter on curve 1
        # is close to the exact one.
        computed[index, (4, 5)] = intersection_info.curve1.evaluate(s_val)
        exact[index, (4, 5)] = point

        # Make sure the point corresponding to the parameter on curve 2
        # is close to the exact one.
        computed[index, (6, 7)] = intersection_info.curve2.evaluate(t_val)
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


def incorrect_values(multipliers, errors, exact):
    observed_error_ulps = np.abs(errors / SPACING(exact))

    failure_rows, failure_cols = np.where(observed_error_ulps > multipliers)
    failures = []
    for failure_row, failure_col in six.moves.zip(failure_rows, failure_cols):
        failures.append('* {:d}, {:d}'.format(failure_row, failure_col))
    failures = '\n'.join(failures)

    zero_misses = np.where((exact == 0.0) & (observed_error_ulps != 0.0))
    observed_error_ulps[zero_misses] = np.nan

    msg = INCORRECT_VALUES_MSG.format(
        multipliers, observed_error_ulps, failures)
    raise AssertionError(msg)


def check_intersect(intersection_info, strategy):
    computed, exact = intersection_values(intersection_info, strategy)
    zero_misses, multipliers = error_multipliers(
        intersection_info, exact.shape, strategy)

    # Make sure our errors fall under the number of "allowed" ULPs.
    allowed_errors = np.abs(multipliers * SPACING(exact))
    errors = np.abs(exact - computed)
    if not np.all(errors <= allowed_errors):
        incorrect_values(multipliers, errors, exact)

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
    tangent_kwargs = TANGENT_OVERRIDES[strategy].get(id_, {})

    if tangent_kwargs == {'success': True}:
        check_intersect(intersection_info, strategy)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            intersection_values(intersection_info, strategy)

        exc_args = exc_info.value.args
        if strategy is ALGEBRAIC:
            assert len(exc_args) == 2
            assert exc_args[0] == NON_SIMPLE_ERR
        else:
            if 'parallel' in tangent_kwargs:
                assert tangent_kwargs == {'parallel': True}
                assert exc_args == SEGMENTS_PARALLEL
            else:
                too_many = tangent_kwargs.get('too_many')
                assert tangent_kwargs == {'too_many': too_many}
                err_msg = TOO_MANY.format(too_many, 4 * too_many)
                assert exc_args == (err_msg,)


def check_coincident(intersection_info, strategy):
    id_ = intersection_info.id_
    coincident_kwargs = COINCIDENT_OVERRIDES[strategy].get(id_, {})

    if coincident_kwargs == {'success': True}:
        check_intersect(intersection_info, strategy)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            intersection_values(intersection_info, strategy)

        exc_args = exc_info.value.args
        if strategy is ALGEBRAIC:
            assert coincident_kwargs == {}
            assert exc_args == COINCIDENT_ERR
        else:
            too_many = coincident_kwargs.get('too_many')
            assert coincident_kwargs == {'too_many': too_many}
            err_msg = TOO_MANY.format(too_many, 4 * too_many)
            assert exc_args == (err_msg,)


@pytest.mark.parametrize(
    'strategy,intersection_info',
    itertools.product(
        (GEOMETRIC, ALGEBRAIC),
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
