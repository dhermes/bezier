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

"""Environment configuration for ``bezier`` runtime dependencies.

Only needed for Windows, to add ``extra-dll`` directory to the DLL search
path so that the ``libbezier`` DLL can be located.
"""

import os

import pkg_resources


# Error messages for ``handle_import_error``.
TEMPLATE = "No module named 'bezier.{}'"  # 3.6, 3.7, 3.8, pypy3
# NOTE: ``os.add_dll_directory()`` was added on Windows in Python 3.8.
OS_ADD_DLL_DIRECTORY = getattr(os, "add_dll_directory", None)


def add_dll_directory(extra_dll_dir):
    """Add a DLL directory.

    This is only expected to be invoked on Windows. For Python versions before
    3.8, this will update the ``%PATH%`` environment variable to include
    ``extra_dll_dir`` and for 3.8 and later, it will invoke
    ``os.add_dll_directory()``.
    """
    if not os.path.isdir(extra_dll_dir):
        return

    if OS_ADD_DLL_DIRECTORY is not None:
        OS_ADD_DLL_DIRECTORY(extra_dll_dir)  # pylint: disable=not-callable
        return

    path = os.environ.get("PATH", "")
    values = [subdir for subdir in path.split(os.pathsep) if subdir]
    values.append(extra_dll_dir)
    os.environ["PATH"] = os.pathsep.join(values)


def modify_path():
    """Add the DLL directory to the module search path.

    This will only modify path if
    * on Windows
    * the ``extra-dll`` directory is in package resources
    """
    if os.name != "nt":
        return

    try:
        extra_dll_dir = pkg_resources.resource_filename("bezier", "extra-dll")
    except ImportError:
        return

    add_dll_directory(extra_dll_dir)


def handle_import_error(caught_exc, name):
    """Allow or re-raise an import error.

    This is to distinguish between expected and unexpected import errors.
    If the module is not found, it simply means the Cython / Fortran speedups
    were not built with the package. If the error message is different, e.g.
    ``... undefined symbol: __curve_intersection_MOD_all_intersections``, then
    the import error **should** be raised.

    Args:
        caught_exc (ImportError): An exception caught when trying to import
            a Cython module.
        name (str): The name of the module. For example, for the module
            ``bezier._curve_speedup``, the name is ``"_curve_speedup"``.

    Raises:
        ImportError: If the error message is different than the basic
            "missing module" error message.
    """
    expected_msg = TEMPLATE.format(name)
    if caught_exc.args == (expected_msg,):
        return

    raise caught_exc


modify_path()
