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

import numpy as np
import pytest

from bezier import _implicitization

import runtime_utils


SPACING = np.spacing  # pylint: disable=no-member
ULPS_ALLOWED = 3.0
# NOTE: We use units of least precision (ULP) as error. These
#       are for the very rare cases where the computed values
#       differ from the actual values by more than 3 ULPs.
CUSTOM_ERRORS = {
    (8, 27): np.asfortranarray([
        [ULPS_ALLOWED, 4.0],
        [ULPS_ALLOWED, 6.0],
    ]),
    (11, 26): np.asfortranarray([
        [12.0, 30.0],
        [ULPS_ALLOWED, ULPS_ALLOWED],
        [ULPS_ALLOWED, ULPS_ALLOWED],
    ]),
    (50, 54): np.asfortranarray([
        [0.0, 165.0],
    ]),
    (51, 54): np.asfortranarray([
        [0.0, 18.0],
    ]),
    (52, 54): np.asfortranarray([
        [0.0, 165.0],
    ]),
    (53, 54): np.asfortranarray([
        [0.0, 18.0],
    ]),
}
CURVES, INTERSECTIONS = runtime_utils.get_intersections_info()


def check_tangent(nodes1, nodes2):
    with pytest.raises(NotImplementedError) as exc_info:
        _implicitization.intersect_curves(nodes1, nodes2)

    assert len(exc_info.value.args) == 2
    assert exc_info.value.args[0] == _implicitization._NON_SIMPLE_ERR


def check_coincident(nodes1, nodes2):
    with pytest.raises(NotImplementedError) as exc_info:
        _implicitization.intersect_curves(nodes1, nodes2)

    assert exc_info.value.args == (_implicitization._COINCIDENT_ERR,)


def check_intersect(id_pair, nodes1, nodes2, info):
    param_vals = _implicitization.intersect_curves(nodes1, nodes2)
    assert param_vals.size > 0

    # NOTE: This assumes the intersections are sorted by s-value.
    s_vals = info['curve1_params']
    num_s, = s_vals.shape
    exact = np.zeros((num_s, 2), order='F')
    exact[:, 0] = s_vals
    exact[:, 1] = info['curve2_params']

    multiplier = CUSTOM_ERRORS.get(id_pair, ULPS_ALLOWED)
    # NOTE: Spacing gives ULP for each value.
    allowed_errors = multiplier * SPACING(exact)

    computed = param_vals[np.argsort(param_vals[:, 0]), :2]
    assert exact.shape == computed.shape
    # NOTE: We assume zeros will be **exactly** correct.
    assert np.all(np.abs(exact - computed) <= allowed_errors)


def check_no_intersect(nodes1, nodes2):
    param_vals = _implicitization.intersect_curves(nodes1, nodes2)
    assert param_vals.size == 0


def id_func(intersection_info):
    return 'curves {:d} and {:d} (ID: {:d})'.format(
        intersection_info['curve1'], intersection_info['curve2'],
        intersection_info['id'])


@pytest.mark.parametrize(
    'intersection_info',
    INTERSECTIONS,
    ids=id_func,
)
def test_intersect(intersection_info):
    # Get info for "curve 1".
    curve_id1 = intersection_info['curve1']
    curve1_info = CURVES[curve_id1 - 1]
    assert curve1_info['id'] == curve_id1
    nodes1 = curve1_info['control_points']

    # Get info for "curve 2".
    curve_id2 = intersection_info['curve2']
    curve2_info = CURVES[curve_id2 - 1]
    assert curve2_info['id'] == curve_id2
    nodes2 = curve2_info['control_points']

    # Actually try to intersect the curves.
    intersection_type = intersection_info['type']
    if intersection_type == 'tangent':
        check_tangent(nodes1, nodes2)
    elif intersection_type == 'coincident':
        check_coincident(nodes1, nodes2)
    elif intersection_type == 'standard':
        id_pair = (curve_id1, curve_id2)
        check_intersect(id_pair, nodes1, nodes2, intersection_info)
    elif intersection_type == 'no-intersection':
        check_no_intersect(nodes1, nodes2)
    else:
        raise ValueError(
            'Unexpected intersection type', intersection_type)
