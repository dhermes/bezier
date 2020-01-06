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

import numpy as np
import pytest
import six

from bezier import _algebraic_intersection
from bezier import _geometric_intersection
from bezier import _py_geometric_intersection
from bezier import _py_intersection_helpers
import bezier.curve
from tests import utils as base_utils
from tests.functional import utils
from tests.functional.utils import CurveIntersectionType


SPACING = np.spacing  # pylint: disable=no-member
GEOMETRIC = bezier.curve.IntersectionStrategy.GEOMETRIC
ALGEBRAIC = bezier.curve.IntersectionStrategy.ALGEBRAIC
_, INTERSECTIONS = utils.curve_intersections_info()
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
ZERO_THRESHOLD = 0.5 ** 52
ZERO_MISS = object()
# NOTE: We use units of least precision (ULP) as error. These
#       are for the very rare cases where the computed values
#       differ from the actual values by more than 3 ULPs.
ULPS_ALLOWED_OVERRIDE = {
    GEOMETRIC: {
        1: {(1, 1): 4},  # Established on CentOS 5 (i686 Docker image)
        12: {(1, 0): 4},  # Established on Ubuntu 16.04
        17: {
            (5, 0): ZERO_MISS,  # Established on Ubuntu 16.04
            (1, 2): 8,  # Established on CentOS 5 (i686 Docker image)
            (1, 3): 4,  # Established on Ubuntu 16.04
        },
        18: {
            (4, 0): 5,  # Established on Ubuntu 16.04
            (0, 1): 8,  # Established on Ubuntu 16.04
            (4, 1): 4,  # Established on Ubuntu 16.04
            (0, 2): 4,  # Established on CentOS 5 (i686 Docker image)
            (3, 3): ZERO_MISS,  # Established on Ubuntu 16.04
            (5, 3): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        22: {
            (0, 0): 8,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 21,  # Established on Ubuntu 16.04
        },
        23: {
            (0, 0): 14,  # Established on Ubuntu 16.04
            (1, 0): 41,  # Established on Ubuntu 16.04
            (0, 1): 16,  # Established on Ubuntu 16.04
            (1, 1): 21,  # Established on Ubuntu 16.04
        },
        25: {(4, 0): ZERO_MISS},  # Established on Ubuntu 16.04
        28: {(2, 0): 4},  # Established on CentOS 5 (i686 Docker image)
        31: {(3, 0): ZERO_MISS},  # Established on CentOS 5 (i686 Docker image)
        37: {(1, 0): 91},  # Established on Ubuntu 16.04
        38: {(1, 0): 1013},  # Established on Ubuntu 16.04
        39: {(1, 0): 91},  # Established on Ubuntu 16.04
        40: {(1, 0): 1013},  # Established on Ubuntu 16.04
        46: {(1, 2): 22},  # Established on Ubuntu 16.04
        49: {
            (0, 0): 13,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 13,  # Established on CentOS 5 (i686 Docker image)
        },
        50: {
            (0, 0): 447,  # Established on Ubuntu 16.04
            (0, 1): 62,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 474,  # Established on Ubuntu 16.04
            (1, 1): 69,  # Established on CentOS 5 (i686 Docker image)
        },
        51: {
            (0, 0): 2259,  # Established on Ubuntu 16.04
            (1, 0): 86670,  # Established on Ubuntu 16.04
        },
        52: {
            (0, 0): 15222,  # Established on Ubuntu 16.04
            (1, 0): 10239,  # Established on Ubuntu 16.04
        },
        53: {
            (0, 0): 285086,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 264170,  # Established on CentOS 5 (i686 Docker image)
        },
    },
    ALGEBRAIC: {
        17: {(5, 0): ZERO_MISS},  # Established on Ubuntu 16.04
        18: {
            (4, 0): 5,  # Established on Ubuntu 16.04
            (4, 1): 4,  # Established on Ubuntu 16.04
            (1, 3): 5,  # Established on CentOS 5 (i686 Docker image)
            (3, 3): ZERO_MISS,  # Established on Ubuntu 16.04
            (5, 3): ZERO_MISS,  # Established on Ubuntu 16.04
        },
        22: {
            (0, 0): 16,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 38,  # Established on CentOS 5 (i686 Docker image)
            (0, 1): 4,  # Established on CentOS 5 (i686 Docker image)
            (1, 1): 4,  # Established on CentOS 5 (i686 Docker image)
        },
        23: {
            (0, 1): 4,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 7,  # Established on Ubuntu 18.04
        },
        25: {(4, 0): ZERO_MISS},  # Established on Ubuntu 16.04
        28: {(2, 0): 4},  # Established on CentOS 5 (i686 Docker image)
        37: {(1, 0): 165},  # Established on Ubuntu 16.04
        38: {(1, 0): 18},  # Established on Ubuntu 16.04
        39: {(1, 0): 165},  # Established on Ubuntu 16.04
        40: {(1, 0): 18},  # Established on Ubuntu 16.04
        49: {
            (0, 0): 27,  # Established on Ubuntu 16.04
            (1, 0): 27,  # Established on Ubuntu 16.04
        },
        50: {
            (0, 0): 528,  # Established on CentOS 5 (i686 Docker image)
            (0, 1): 75,  # Established on Ubuntu 18.04
            (1, 0): 551,  # Established on CentOS 5 (i686 Docker image)
            (1, 1): 81,  # Established on Ubuntu 18.04
        },
        51: {
            (0, 0): 3274,  # Established on CentOS 5 (i686 Docker image)
            (1, 0): 120841,  # Established on CentOS 5 (i686 Docker image)
        },
        52: {
            (0, 0): 13677,  # Established on Ubuntu 16.04
            (1, 0): 9260,  # Established on Ubuntu 16.04
        },
    },
}
NON_SIMPLE_ERR = _algebraic_intersection._NON_SIMPLE_ERR
TOO_MANY = _py_geometric_intersection._TOO_MANY_TEMPLATE
BAD_MULTIPLICITY = (_py_intersection_helpers.NEWTON_NO_CONVERGE,)
COINCIDENT_ERR = (_algebraic_intersection._COINCIDENT_ERR,)
TANGENT_OVERRIDES = {
    GEOMETRIC: {
        4: {"success": True},
        11: {"success": True},
        14: {"success": True},
        19: {"success": True},
        24: {"success": True},
        31: {"success": True},
        41: {"success": True},
        42: {"bad_multiplicity": True},
        43: {"success": True},
        44: {"success": True},
        45: {"too_many": 74},
        46: {"success": True},
        47: {"success": True},
        50: {"success": True},
    },
    ALGEBRAIC: {},
}
COINCIDENT_OVERRIDES = {
    GEOMETRIC: {
        20: {"success": True},
        33: {"success": True},
        34: {"success": True},
        35: {"success": True},
    },
    ALGEBRAIC: {53: {}},
}
INCORRECT_COUNT = {GEOMETRIC: (), ALGEBRAIC: ()}
if base_utils.IS_PYPY:
    INCORRECT_COUNT[ALGEBRAIC] += (10,)


def get_sorted_intersections(intersection_info, strategy):
    nodes1 = intersection_info.nodes1
    nodes2 = intersection_info.nodes2
    if strategy is GEOMETRIC:
        intersections, _ = _geometric_intersection.all_intersections(
            nodes1, nodes2
        )
    else:
        intersections, _ = _algebraic_intersection.all_intersections(
            nodes1, nodes2
        )
    # Make we have the right number of intersections.
    if intersections.shape != (2, intersection_info.num_params):
        raise utils.IncorrectCount(
            "Received wrong number of intersections",
            intersections.shape,
            "Expected",
            intersection_info.num_params,
            intersection_info.test_id,
        )

    # Sort the intersections by s-value.
    sorted_inds = np.argsort(intersections[0, :])
    return intersections[:, sorted_inds]


def intersection_values(intersection_info, strategy):
    # Actually intersect the curves.
    intersections = get_sorted_intersections(intersection_info, strategy)
    # Put the computed and expected values into matrices to be compared.
    s_vals, t_vals, intersection_pts = intersection_info.params
    computed = np.zeros((6, intersection_info.num_params), order="F")
    exact = np.zeros(computed.shape, order="F")
    info = six.moves.zip(intersections.T, s_vals, t_vals, intersection_pts.T)
    for index, (intersection, s_val, t_val, point) in enumerate(info):
        computed[(0, 1), index] = intersection
        exact[(0, 1), index] = s_val, t_val
        # Make sure the point corresponding to the parameter on curve 1
        # is close to the exact one.
        computed[(2, 3), index] = intersection_info.curve1.evaluate(s_val)[
            :, 0
        ]
        exact[(2, 3), index] = point
        # Make sure the point corresponding to the parameter on curve 2
        # is close to the exact one.
        computed[(4, 5), index] = intersection_info.curve2.evaluate(t_val)[
            :, 0
        ]
        exact[(4, 5), index] = point
    return computed, exact


def error_multipliers(intersection_info, shape, strategy):
    zero_misses = []
    multipliers = ULPS_ALLOWED * np.ones(shape, order="F")
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
    # pylint: disable=unbalanced-tuple-unpacking
    failure_rows, failure_cols = np.where(observed_error_ulps > multipliers)
    # pylint: enable=unbalanced-tuple-unpacking
    failures = []
    for failure_row, failure_col in six.moves.zip(failure_rows, failure_cols):
        failures.append("* {:d}, {:d}".format(failure_row, failure_col))
    failures = "\n".join(failures)
    zero_misses = np.where((exact == 0.0) & (observed_error_ulps != 0.0))
    observed_error_ulps[zero_misses] = np.nan
    msg = INCORRECT_VALUES_MSG.format(
        multipliers, observed_error_ulps, failures
    )
    raise AssertionError(msg)


def check_intersect(intersection_info, strategy):
    computed, exact = intersection_values(intersection_info, strategy)
    zero_misses, multipliers = error_multipliers(
        intersection_info, exact.shape, strategy
    )
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
    if tangent_kwargs == {"success": True}:
        check_intersect(intersection_info, strategy)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            intersection_values(intersection_info, strategy)
        exc_args = exc_info.value.args
        if strategy is ALGEBRAIC:
            assert len(exc_args) == 2
            assert exc_args[0] == NON_SIMPLE_ERR
        else:
            if "bad_multiplicity" in tangent_kwargs:
                assert tangent_kwargs == {"bad_multiplicity": True}
                assert exc_args == BAD_MULTIPLICITY
            else:
                too_many = tangent_kwargs.get("too_many")
                assert tangent_kwargs == {"too_many": too_many}
                err_msg = TOO_MANY.format(too_many)
                assert exc_args == (err_msg,)


def check_coincident(intersection_info, strategy):
    id_ = intersection_info.id_
    coincident_kwargs = COINCIDENT_OVERRIDES[strategy].get(id_, {})
    if coincident_kwargs == {"success": True}:
        check_intersect(intersection_info, strategy)
    else:
        with pytest.raises(NotImplementedError) as exc_info:
            intersection_values(intersection_info, strategy)
        exc_args = exc_info.value.args
        if strategy is ALGEBRAIC:
            assert coincident_kwargs == {}
            assert exc_args == COINCIDENT_ERR
        else:
            too_many = coincident_kwargs.get("too_many")
            assert coincident_kwargs == {"too_many": too_many}
            err_msg = TOO_MANY.format(too_many)
            assert exc_args == (err_msg,)


@pytest.mark.parametrize(
    "strategy,intersection_info",
    itertools.product((GEOMETRIC, ALGEBRAIC), INTERSECTIONS),
    ids=utils.id_func,
)
def test_intersect(strategy, intersection_info):
    id_ = intersection_info.id_
    # Actually try to intersect the curves.
    intersection_type = intersection_info.type_
    if id_ in INCORRECT_COUNT[strategy]:
        with pytest.raises(utils.IncorrectCount):
            check_intersect(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.tangent:
        check_tangent(intersection_info, strategy)
    elif (
        intersection_type == CurveIntersectionType.coincident
        or id_ in COINCIDENT_OVERRIDES[strategy]
    ):
        check_coincident(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.standard:
        check_intersect(intersection_info, strategy)
    elif intersection_type == CurveIntersectionType.no_intersection:
        check_no_intersect(intersection_info, strategy)
    else:
        raise ValueError("Unexpected intersection type", intersection_type)
