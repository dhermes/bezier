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

Only needed for Windows, to add extra directory to the DLL search
path so that the ``libbezier`` DLL can be located.
"""

import os


# Error messages for ``handle_import_error``.
TEMPLATE = "No module named 'bezier.{}'"  # 3.10, 3.11, 3.12, pypy3
EXTRA_DLL_ENV = "BEZIER_EXTRA_DLL"
"""Environment variable used to add extra directory to DLL search path.

This is intended to be used in tests and when building from source.
"""


def add_dll_directory():
    """Add a DLL directory.

    This is only expected to be invoked on Windows. It will invoke
    ``os.add_dll_directory()``, which was added in Python 3.8.
    """
    if os.name != "nt":
        return

    bezier_extra_dll = os.environ.get(EXTRA_DLL_ENV)
    if bezier_extra_dll is None:
        return

    extra_dll_dirs = bezier_extra_dll.split(os.pathsep)
    for extra_dll_dir in extra_dll_dirs:
        if not os.path.isdir(extra_dll_dir):
            continue

        os.add_dll_directory(extra_dll_dir)


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


add_dll_directory()
