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

"""Helpers for intersecting B |eacute| zier curves via geometric methods."""

import atexit

from bezier import _py_geometric_intersection

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
Linearization = _py_geometric_intersection.Linearization
BoxIntersectionType = _py_geometric_intersection.BoxIntersectionType
segment_intersection = _py_geometric_intersection.segment_intersection

if _speedup is None:  # pragma: NO COVER
    bbox_intersect = _py_geometric_intersection.bbox_intersect
    all_intersections = _py_geometric_intersection.all_intersections
else:
    bbox_intersect = _speedup.bbox_intersect
    all_intersections = _speedup.curve_intersections
    atexit.register(_speedup.free_curve_intersections_workspace)
# pylint: enable=invalid-name
