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

"""Helpers for ``setup.py`` specific to OS X."""

import sys

import setup_helpers


MAC_OS_X = "darwin"


def is_osx_gfortran(f90_compiler):
    """Checks if the current build is ``gfortran`` on OS X.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.

    Returns:
        bool: Only :data:`True` if

        * Current OS is OS X (checked via ``sys.platform``).
        * ``f90_compiler`` corresponds to ``gfortran``.
    """
    # NOTE: NumPy may not be installed, but we don't want **this** module to
    #       cause an import failure.
    from numpy.distutils.fcompiler import gnu

    # Only Mac OS X.
    if sys.platform != MAC_OS_X:
        return False

    # Only ``gfortran``.
    if not isinstance(f90_compiler, gnu.Gnu95FCompiler):
        return False

    return True


def patch_f90_compiler(f90_compiler):
    """Patch up ``f90_compiler.library_dirs``.

    On Mac OS X, a Homebrew installed ``gfortran`` needs some help. The
    ``numpy.distutils`` "default" constructor for ``Gnu95FCompiler`` only has
    a single library search path, but there are many library paths included in
    the full ``gcc`` install.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.
    """
    if not is_osx_gfortran(f90_compiler):
        return

    library_dirs = f90_compiler.library_dirs
    # ``library_dirs`` is a list (i.e. mutable), so we can update in place.
    library_dirs[:] = setup_helpers.gfortran_search_path(library_dirs)
