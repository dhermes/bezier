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

import numpy as np
import pytest

from bezier import _implicitization

import runtime_utils
from runtime_utils import CurveIntersectionType


SPACING = np.spacing  # pylint: disable=no-member
ULPS_ALLOWED = 3.0
# NOTE: We use units of least precision (ULP) as error. These
#       are for the very rare cases where the computed values
#       differ from the actual values by more than 3 ULPs.
CUSTOM_ERRORS = {
    22: np.asfortranarray([
        [12.0, 30.0],
        [ULPS_ALLOWED, ULPS_ALLOWED],
        [ULPS_ALLOWED, ULPS_ALLOWED],
    ]),
    23: np.asfortranarray([
        [ULPS_ALLOWED, 4.0],
        [ULPS_ALLOWED, 6.0],
    ]),
    37: np.asfortranarray([
        [0.0, 165.0],
    ]),
    38: np.asfortranarray([
        [0.0, 18.0],
    ]),
    39: np.asfortranarray([
        [0.0, 165.0],
    ]),
    40: np.asfortranarray([
        [0.0, 18.0],
    ]),
}
_, INTERSECTIONS = runtime_utils.curve_intersections_info()


def check_tangent(nodes1, nodes2):
    with pytest.raises(NotImplementedError) as exc_info:
        _implicitization.intersect_curves(nodes1, nodes2)

    assert len(exc_info.value.args) == 2
    assert exc_info.value.args[0] == _implicitization._NON_SIMPLE_ERR


def check_coincident(nodes1, nodes2):
    with pytest.raises(NotImplementedError) as exc_info:
        _implicitization.intersect_curves(nodes1, nodes2)

    assert exc_info.value.args == (_implicitization._COINCIDENT_ERR,)


def check_intersect(nodes1, nodes2, intersection_info):
    param_vals = _implicitization.intersect_curves(nodes1, nodes2)
    assert param_vals.size > 0

    # NOTE: This assumes the intersections are sorted by s-value.
    exact = np.zeros((intersection_info.num_params, 2), order='F')
    exact[:, 0] = intersection_info.curve1_params
    exact[:, 1] = intersection_info.curve2_params

    multiplier = CUSTOM_ERRORS.get(intersection_info.id_, ULPS_ALLOWED)
    # NOTE: Spacing gives ULP for each value.
    allowed_errors = multiplier * SPACING(exact)

    computed = param_vals[np.argsort(param_vals[:, 0]), :]
    assert exact.shape == computed.shape
    # NOTE: We assume zeros will be **exactly** correct.
    assert np.all(np.abs(exact - computed) <= allowed_errors)


def check_no_intersect(nodes1, nodes2):
    param_vals = _implicitization.intersect_curves(nodes1, nodes2)
    assert param_vals.size == 0


@pytest.mark.parametrize(
    'intersection_info',
    INTERSECTIONS,
    ids=operator.attrgetter('tests_id'),
)
def test_intersect(intersection_info):
    # Get the control points for the curves.
    nodes1 = intersection_info.curve1.control_points
    nodes2 = intersection_info.curve2.control_points

    # Actually try to intersect the curves.
    intersection_type = intersection_info.type_
    if intersection_type == CurveIntersectionType.tangent:
        check_tangent(nodes1, nodes2)
    elif intersection_type == CurveIntersectionType.coincident:
        check_coincident(nodes1, nodes2)
    elif intersection_type == CurveIntersectionType.standard:
        check_intersect(nodes1, nodes2, intersection_info)
    elif intersection_type == CurveIntersectionType.no_intersection:
        check_no_intersect(nodes1, nodes2)
    else:
        raise ValueError(
            'Unexpected intersection type', intersection_type)
