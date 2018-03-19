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

.. autoclass:: UnsupportedDegree
   :members:
"""

import os

import pkg_resources

# NOTE: ``__config__`` **must** be the first import because it (may)
#       modify the search path used to locate shared libraries.
from bezier import __config__
from bezier._helpers import UnsupportedDegree
from bezier.curve import Curve
from bezier.curved_polygon import CurvedPolygon
from bezier.surface import Surface

try:
    import bezier._speedup  # noqa: F401

    _HAS_SPEEDUP = True
except ImportError as exc:  # pragma: NO COVER
    __config__.handle_import_error(exc, '_speedup')
    _HAS_SPEEDUP = False
# NOTE: The ``__version__`` and ``__author__`` are hard-coded here, rather
#       than using ``pkg_resources.get_distribution('bezier').version``
#       and related. This is **entirely** to accomodate builds where
#       ``bezier`` is imported from source (and not installed).
__author__ = 'Danny Hermes'
__version__ = '0.8.0'
"""str: The current version of :mod:`bezier`."""
__all__ = [
    '__author__',
    '__version__',
    'Curve',
    'CurvedPolygon',
    'get_dll',
    'get_include',
    'get_lib',
    'Surface',
    'UnsupportedDegree',
]


def get_include():
    """Get the directory with ``.h`` header files.

    Extension modules (and Cython modules) that need to compile against
    ``libbezier`` should use this function to locate the appropriate
    include directory.

    For more information, see :doc:`../native-libraries`.

    Returns:
        str: ``include`` directory that contains header files for the
        ``libbezier`` Fortran library.
    """
    return pkg_resources.resource_filename('bezier', 'include')


def get_lib():
    """Get the directory with ``.a`` / ``.lib`` static libraries.

    Extension modules (and Cython modules) that need to compile against
    ``libbezier`` should use this function to locate the appropriate lib
    directory.

    For more information, see :doc:`../native-libraries`.

    Returns:
        str: ``lib`` directory that contains static libraries for the
        ``libbezier`` Fortran library.
    """
    return pkg_resources.resource_filename('bezier', 'lib')


def get_dll():
    """Get the directory with the Windows shared library.

    Extension modules (and Cython modules) that need to compile against
    ``libbezier`` should use this function to locate the appropriate
    Windows shared library or libraries (DLLs).

    For more information, see :doc:`../native-libraries`.

    Returns:
        str: ``extra-dll`` directory that contains the Windows shared library
        for the ``libbezier`` Fortran library.

    Raises:
        OSError: If this function is used anywhere other than Windows.
    """
    if os.name == 'nt':
        return pkg_resources.resource_filename('bezier', 'extra-dll')

    else:
        raise OSError('This function should only be used on Windows.')
