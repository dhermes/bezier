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

"""Helper methods for :mod:`bezier.triangle`.

The functions provided by this module have a Cython speedup with the
exact same interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built.
"""

from bezier.hazmat import triangle_helpers

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    de_casteljau_one_round = triangle_helpers.de_casteljau_one_round
    specialize_triangle = triangle_helpers.specialize_triangle
    subdivide_nodes = triangle_helpers.subdivide_nodes
    jacobian_both = triangle_helpers.jacobian_both
    jacobian_det = triangle_helpers.jacobian_det
    evaluate_barycentric = triangle_helpers.evaluate_barycentric
    evaluate_barycentric_multi = triangle_helpers.evaluate_barycentric_multi
    evaluate_cartesian_multi = triangle_helpers.evaluate_cartesian_multi
    compute_edge_nodes = triangle_helpers.compute_edge_nodes
    compute_area = triangle_helpers.compute_area
else:
    de_casteljau_one_round = _speedup.de_casteljau_one_round
    specialize_triangle = _speedup.specialize_triangle
    subdivide_nodes = _speedup.subdivide_nodes_triangle
    jacobian_both = _speedup.jacobian_both
    jacobian_det = _speedup.jacobian_det
    evaluate_barycentric = _speedup.evaluate_barycentric
    evaluate_barycentric_multi = _speedup.evaluate_barycentric_multi
    evaluate_cartesian_multi = _speedup.evaluate_cartesian_multi
    compute_edge_nodes = _speedup.compute_edge_nodes
    compute_area = _speedup.compute_area
# pylint: enable=invalid-name
