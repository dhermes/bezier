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

.. testsetup:: show-headers, show-lib, show-pxd

   import os
   import textwrap

   import bezier

   def sort_key(name):
       return name.lower().lstrip('_')

   def tree(directory, suffix=None):
       names = sorted(os.listdir(directory), key=sort_key)
       parts = []
       for name in names:
           path = os.path.join(directory, name)
           if os.path.isdir(path):
               sub_part = tree(path, suffix=suffix)
               if sub_part is not None:
                   parts.append(name + os.path.sep)
                   parts.append(textwrap.indent(sub_part, '  '))
           else:
               if suffix is None or name.endswith(suffix):
                   parts.append(name)

       if parts:
           return '\n'.join(parts)
       else:
           return None

   def print_tree(directory, suffix=None):
       print(os.path.basename(directory) + os.path.sep)
       full_tree = tree(directory, suffix=suffix)
       print(textwrap.indent(full_tree, '  '))
"""

import pkg_resources

from bezier.curve import Curve
from bezier.curved_polygon import CurvedPolygon
from bezier.surface import Surface
try:
    import bezier._curve_intersection_speedup  # noqa: F401
    import bezier._curve_speedup  # noqa: F401
    import bezier._helpers_speedup  # noqa: F401
    import bezier._surface_speedup  # noqa: F401
    _HAS_SPEEDUP = True
except ImportError:  # pragma: NO COVER
    _HAS_SPEEDUP = False


__all__ = [
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

    For example:

    .. testsetup:: setup-extension

       import bezier

    .. doctest:: setup-extension

       >>> import setuptools
       >>>
       >>> extension = setuptools.Extension(
       ...     'wrapper',
       ...     ['wrapper.c'],
       ...     include_dirs=[
       ...         bezier.get_include(),
       ...     ],
       ...     libraries=['bezier'],
       ...     library_dirs=[
       ...         bezier.get_lib(),
       ...     ],
       ... )
       >>> extension
       <setuptools.extension.Extension('wrapper') at 0x...>

    The headers are in the ``bezier`` subdirectory and there is a
    catch-all ``bezier.h`` that just includes all of the headers:

    .. doctest:: show-headers

       >>> include_directory = bezier.get_include()
       >>> include_directory
       '.../site-packages/bezier/include'
       >>> print_tree(include_directory)
       include/
         bezier/
           curve.h
           curve_intersection.h
           helpers.h
           surface.h
         bezier.h

    In addition to the header files, several ``cimport``-able ``.pxd``
    Cython declaration files are provided:

    .. doctest:: show-pxd

       >>> include_directory = bezier.get_include()
       >>> bezier_directory = os.path.dirname(include_directory)
       >>> bezier_directory
       '.../site-packages/bezier'
       >>> print_tree(bezier_directory, suffix='.pxd')
       bezier/
         _curve.pxd
         _curve_intersection.pxd
         _helpers.pxd
         _surface.pxd

    For example, ``cimport bezier._curve`` will provide all the functions
    in ``bezier/curve.h``.

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

    See :func:`~bezier.get_include` for an example on how / when to use
    this function.

    We expect a single static library (a ``.lib`` file on Windows and a
    ``.a`` file elsewhere):

    .. doctest:: show-lib

       >>> lib_directory = bezier.get_lib()
       >>> lib_directory
       '.../site-packages/bezier/lib'
       >>> print_tree(lib_directory)
       lib/
         libbezier.a

    Returns:
        str: ``lib`` directory that contains static libraries for the
        ``libbezier`` Fortran library.
    """
    return pkg_resources.resource_filename('bezier', 'lib')
