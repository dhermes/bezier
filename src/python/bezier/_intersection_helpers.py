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

"""Helper methods for intersecting B |eacute| zier shapes.

The functions provided by this module have a Cython speedup with the
exact same interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built.
"""

from bezier.hazmat import intersection_helpers

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    newton_refine = intersection_helpers.newton_refine
else:
    newton_refine = _speedup.newton_refine_curve_intersect
# pylint: enable=invalid-name
