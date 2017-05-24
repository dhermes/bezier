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
import six

from bezier import _implicitization

import candidate_curves


SPACING = np.spacing  # pylint: disable=no-member
# NOTE: This is much too large (for now) but will be resolved
#       in subsequent iterations.
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
        [10.0, 30.0],
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
TANGENT_INTERSECTIONS = (
    (1, 6),
    (10, 23),
    (14, 15),
    (28, 29),
    (38, 39),
    (1, 18),
)
COINCIDENT_INTERSECTIONS = (
    (1, 24),
    (42, 43),
    (44, 45),
    (46, 47),
)


def test_all():
    all_intersect = six.iteritems(candidate_curves.INTERSECTION_INFO)
    for (curve_id1, curve_id2), info in all_intersect:
        id_pair = (curve_id1, curve_id2)
        curve1 = candidate_curves.CURVES[curve_id1]
        nodes1 = curve1._nodes
        curve2 = candidate_curves.CURVES[curve_id2]
        nodes2 = curve2._nodes
        if id_pair in TANGENT_INTERSECTIONS:
            with pytest.raises(NotImplementedError) as exc_info:
                _implicitization.intersect_curves(nodes1, nodes2)

            assert len(exc_info.value.args) == 2
            assert exc_info.value.args[0] == _implicitization._NON_SIMPLE_ERR
            continue
        elif id_pair in COINCIDENT_INTERSECTIONS:
            with pytest.raises(NotImplementedError) as exc_info:
                _implicitization.intersect_curves(nodes1, nodes2)

            assert exc_info.value.args == (_implicitization._COINCIDENT_ERR,)
            continue

        param_vals = _implicitization.intersect_curves(nodes1, nodes2)
        if param_vals.size == 0:
            assert info.size == 0
        else:
            # NOTE: This assumes ``info`` is sorted by s-value.
            exact = info[:, :2]
            multiplier = CUSTOM_ERRORS.get(id_pair, ULPS_ALLOWED)
            # NOTE: Spacing gives ULP for each value.
            allowed_errors = multiplier * SPACING(exact)

            computed = param_vals[np.argsort(param_vals[:, 0]), :2]
            assert exact.shape == computed.shape
            # NOTE: We assume zeros will be **exactly** correct.
            assert np.all(np.abs(exact - computed) <= allowed_errors)
