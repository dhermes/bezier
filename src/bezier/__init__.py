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

r"""Helper for B |eacute| zier Curves, Triangles, and Higher Order Objects.

Intended to perform basic operations on B |eacute| zier objects such
as intersections, length/area/etc. computations, subdivision,
implicitization and other relevant information.

Plotting utilities are also provided.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import pkg_resources

# NOTE: ``__config__`` **must** be the first import because it (may)
#       modify the search path used to locate shared libraries.
from bezier import __config__  # noqa: F401
from bezier.curve import Curve
from bezier.curved_polygon import CurvedPolygon
from bezier.surface import Surface
try:
    import bezier._helpers_speedup  # noqa: F401
    _HAS_HELPERS_SPEEDUP = True
except ImportError:  # pragma: NO COVER
    _HAS_HELPERS_SPEEDUP = False
try:
    import bezier._curve_speedup  # noqa: F401
    _HAS_CURVE_SPEEDUP = True
except ImportError:  # pragma: NO COVER
    _HAS_CURVE_SPEEDUP = False
try:
    import bezier._surface_speedup  # noqa: F401
    _HAS_SURFACE_SPEEDUP = True
except ImportError:  # pragma: NO COVER
    _HAS_SURFACE_SPEEDUP = False
try:
    import bezier._curve_intersection_speedup  # noqa: F401
    _HAS_CURVE_INTERSECTION_SPEEDUP = True
except ImportError:  # pragma: NO COVER
    _HAS_CURVE_INTERSECTION_SPEEDUP = False


__version__ = pkg_resources.get_distribution('bezier').version
"""str: The current version of :mod:`bezier`."""
__all__ = [
    '__version__',
    'Curve',
    'CurvedPolygon',
    'get_include',
    'get_lib',
    'Surface',
]


def get_include():
    """Get the directory with ``.h`` header files.

    Extension modules (and Cython modules) that need to compile against
    ``bezier`` should use this function to locate the appropriate include
    directory.

    For more information, see :doc:`../native-libraries`.

    Returns:
        str: ``include`` directory that contains header files for the
        ``libbezier`` Fortran library.
    """
    return pkg_resources.resource_filename('bezier', 'include')


def get_lib():
    """Get the directory with ``.a`` / ``.lib`` static libraries.

    Extension modules (and Cython modules) that need to compile against
    the ``libbezier`` static library should use this function to locate the
    appropriate lib directory.

    For more information, see :doc:`../native-libraries`.

    Returns:
        str: ``lib`` directory that contains static libraries for the
        ``libbezier`` Fortran library.
    """
    return pkg_resources.resource_filename('bezier', 'lib')
