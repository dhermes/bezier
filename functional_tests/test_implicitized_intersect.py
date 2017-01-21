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
import six

from bezier import _implicitization

import candidate_curves


# NOTE: This is much too large (for now) but will be resolved
#       in subsequent iterations.
EPS = 0.5**45
IGNORED = (
    (1, 13),  # Degree pair (1-4) unsupported.
    (1, 24),  # Coincident
    (15, 25),  # Spurious root
    (42, 43),  # Coincident
    (44, 45),  # Coincident
    (46, 47),  # Coincident
)


def test_all():
    all_intersect = six.iteritems(candidate_curves.INTERSECTION_INFO)
    for (curve_id1, curve_id2), info in all_intersect:
        if (curve_id1, curve_id2) in IGNORED:
            continue
        curve1 = candidate_curves.CURVES[curve_id1]
        nodes1 = curve1._nodes
        curve2 = candidate_curves.CURVES[curve_id2]
        nodes2 = curve2._nodes
        param_vals = _implicitization.intersect_curves(nodes1, nodes2)
        if param_vals.size == 0:
            assert info.size == 0
        else:
            if param_vals[0, 0] == -np.inf:
                col_id = 1
            else:
                col_id = 0

            exact = np.sort(info[:, col_id])
            computed = np.sort(param_vals[:, col_id])
            assert np.allclose(exact, computed,
                               atol=0.0, rtol=EPS)
