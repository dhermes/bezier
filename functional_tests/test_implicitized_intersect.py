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

import numpy as np
import pytest
import six

from bezier import _implicitization

import candidate_curves


# NOTE: This is much too large (for now) but will be resolved
#       in subsequent iterations.
EPS = 0.5**48
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
        curve1 = candidate_curves.CURVES[curve_id1]
        nodes1 = curve1._nodes
        curve2 = candidate_curves.CURVES[curve_id2]
        nodes2 = curve2._nodes
        if (curve_id1, curve_id2) in TANGENT_INTERSECTIONS:
            with pytest.raises(NotImplementedError) as exc_info:
                _implicitization.intersect_curves(nodes1, nodes2)

            assert len(exc_info.value.args) == 4
            assert exc_info.value.args[0] == _implicitization._NON_SIMPLE_ERR
            continue
        elif (curve_id1, curve_id2) in COINCIDENT_INTERSECTIONS:
            with pytest.raises(NotImplementedError) as exc_info:
                _implicitization.intersect_curves(nodes1, nodes2)

            assert exc_info.value.args == (_implicitization._COINCIDENT_ERR,)
            continue

        param_vals = _implicitization.intersect_curves(nodes1, nodes2)
        if param_vals.size == 0:
            assert info.size == 0
        else:
            exact = info[np.argsort(info[:, 0]), :2]
            computed = param_vals[np.argsort(param_vals[:, 0]), :2]
            zero_elts = np.where(exact == 0.0)
            assert np.all(np.abs(computed[zero_elts]) < EPS)
            computed[zero_elts] = 0.0  # So we can still use atol=0.0.
            assert np.allclose(exact, computed,
                               atol=0.0, rtol=EPS)
