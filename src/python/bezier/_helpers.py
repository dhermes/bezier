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

"""Generic geometry and floating point helpers.

The functions provided by this module have a Cython speedup with the
exact same interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built.
"""

from bezier.hazmat import helpers as _py_helpers

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    vector_close = _py_helpers.vector_close
    in_interval = _py_helpers.in_interval
    bbox = _py_helpers.bbox
    contains_nd = _py_helpers.contains_nd
    cross_product = _py_helpers.cross_product
    wiggle_interval = _py_helpers.wiggle_interval
    simple_convex_hull = _py_helpers.simple_convex_hull
    polygon_collide = _py_helpers.polygon_collide
else:
    vector_close = _speedup.vector_close
    in_interval = _speedup.in_interval
    bbox = _speedup.bbox
    contains_nd = _speedup.contains_nd
    cross_product = _speedup.cross_product
    wiggle_interval = _speedup.wiggle_interval
    simple_convex_hull = _speedup.simple_convex_hull
    polygon_collide = _speedup.polygon_collide
# pylint: enable=invalid-name
