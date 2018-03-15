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

Only needed for Windows, to add ``extra-dll`` directory to the search
``%PATH%`` so that the ``libbezier`` DLL can be located.
"""

import os

import pkg_resources

# Error messages for ``handle_import_error``.
TEMPLATES = (
    'No module named \'bezier.{}\'',  # 3.5, 3.6, pypy3
    'No module named {}',  # 2.7
    'No module named bezier.{}',  # pypy2
)


def modify_path():
    """Modify the module search path."""
    # Only modify path on Windows.
    if os.name != 'nt':
        return

    path = os.environ.get('PATH')
    if path is None:
        return

    try:
        extra_dll_dir = pkg_resources.resource_filename('bezier', 'extra-dll')
        if os.path.isdir(extra_dll_dir):
            os.environ['PATH'] = path + os.pathsep + extra_dll_dir
    except ImportError:
        pass


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
            ``bezier._curve_speedup``, the name is ``'_curve_speedup'``.

    Raises:
        ImportError: If the error message is different than the basic
            "missing module" error message.
    """
    for template in TEMPLATES:
        expected_msg = template.format(name)
        if caught_exc.args == (expected_msg,):
            return

    raise caught_exc


modify_path()
