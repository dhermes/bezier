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

"""Helper methods for B |eacute| zier curves.

The functions provided by this module have a Cython speedup with the
exact same interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built.
"""

from bezier.hazmat import curve_helpers

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    subdivide_nodes = curve_helpers.subdivide_nodes
    evaluate_multi = curve_helpers.evaluate_multi
    evaluate_multi_barycentric = curve_helpers.evaluate_multi_barycentric
    compute_length = curve_helpers.compute_length
    elevate_nodes = curve_helpers.elevate_nodes
    specialize_curve = curve_helpers.specialize_curve
    evaluate_hodograph = curve_helpers.evaluate_hodograph
    get_curvature = curve_helpers.get_curvature
    newton_refine = curve_helpers.newton_refine
    locate_point = curve_helpers.locate_point
    reduce_pseudo_inverse = curve_helpers.reduce_pseudo_inverse
    full_reduce = curve_helpers.full_reduce
else:
    subdivide_nodes = _speedup.subdivide_nodes_curve
    evaluate_multi = _speedup.evaluate_multi
    evaluate_multi_barycentric = _speedup.evaluate_multi_barycentric
    compute_length = _speedup.compute_length
    elevate_nodes = _speedup.elevate_nodes
    specialize_curve = _speedup.specialize_curve
    evaluate_hodograph = _speedup.evaluate_hodograph
    get_curvature = _speedup.get_curvature
    newton_refine = _speedup.newton_refine_curve
    locate_point = _speedup.locate_point_curve
    reduce_pseudo_inverse = _speedup.reduce_pseudo_inverse
    full_reduce = _speedup.full_reduce
# pylint: enable=invalid-name
