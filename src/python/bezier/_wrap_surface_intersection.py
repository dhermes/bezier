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

import atexit

from bezier import _py_surface_intersection

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
algebraic_intersect = _py_surface_intersection.algebraic_intersect

if _speedup is None:  # pragma: NO COVER
    newton_refine = _py_surface_intersection.newton_refine
    locate_point = _py_surface_intersection.locate_point
    geometric_intersect = _py_surface_intersection.geometric_intersect
else:
    newton_refine = _speedup.newton_refine_surface
    locate_point = _speedup.locate_point_surface
    geometric_intersect = _speedup.surface_intersections
    atexit.register(_speedup.free_surface_intersections_workspace)
# pylint: enable=invalid-name
