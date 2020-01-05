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

"""Private helper methods for :mod:`bezier.surface`."""

from bezier import _py_surface_helpers

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
quadratic_jacobian_polynomial = (
    _py_surface_helpers.quadratic_jacobian_polynomial
)
polynomial_sign = _py_surface_helpers.polynomial_sign
cubic_jacobian_polynomial = _py_surface_helpers.cubic_jacobian_polynomial
polynomial_sign = _py_surface_helpers.polynomial_sign

if _speedup is None:  # pragma: NO COVER
    de_casteljau_one_round = _py_surface_helpers.de_casteljau_one_round
    specialize_surface = _py_surface_helpers.specialize_surface
    subdivide_nodes = _py_surface_helpers.subdivide_nodes
    jacobian_both = _py_surface_helpers.jacobian_both
    jacobian_det = _py_surface_helpers.jacobian_det
    evaluate_barycentric = _py_surface_helpers.evaluate_barycentric
    evaluate_barycentric_multi = _py_surface_helpers.evaluate_barycentric_multi
    evaluate_cartesian_multi = _py_surface_helpers.evaluate_cartesian_multi
    compute_edge_nodes = _py_surface_helpers.compute_edge_nodes
    compute_area = _py_surface_helpers.compute_area
else:
    de_casteljau_one_round = _speedup.de_casteljau_one_round
    specialize_surface = _speedup.specialize_surface
    subdivide_nodes = _speedup.subdivide_nodes_surface
    jacobian_both = _speedup.jacobian_both
    jacobian_det = _speedup.jacobian_det
    evaluate_barycentric = _speedup.evaluate_barycentric
    evaluate_barycentric_multi = _speedup.evaluate_barycentric_multi
    evaluate_cartesian_multi = _speedup.evaluate_cartesian_multi
    compute_edge_nodes = _speedup.compute_edge_nodes
    compute_area = _speedup.compute_area
# pylint: enable=invalid-name
