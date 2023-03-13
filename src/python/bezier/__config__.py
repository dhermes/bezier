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


# Error messages for ``handle_import_error``.
TEMPLATE = "No module named 'bezier.{}'"  # 3.8, 3.9, 3.10, 3.11, pypy3


def add_dll_directory(extra_dll_dir):
    """Add a DLL directory.

    This is only expected to be invoked on Windows. It will invoke
    ``os.add_dll_directory()``, which was added in Python 3.8.

    Args:
        extra_dll_dir (str): The path to a directory ``extra-dll``.
    """
    if not os.path.isdir(extra_dll_dir):
        return

    os.add_dll_directory(extra_dll_dir)


def _is_extra_dll(path):
    """Determine if a package path is the extra DLL on Windows.

    Args:
        path (importlib.metadata.PackagePath): A package path.

    Returns:
        bool: Indicating if this is the extra DLL on Windows.
    """
    return "extra-dll" in path.parts and path.name.endswith(".dll")


def _get_extra_dll_dir(bezier_files):
    """Determine if a package path is the extra DLL on Windows.

    Args:
        bezier_files (List[importlib.metadata.PackagePath]): List of package
            paths.

    Returns:
        Optional[str]: The path of the matching ``extra-dll`` directory, or
        :data:`None` if no match can be found.
    """
    for path in bezier_files:
        if not _is_extra_dll(path):
            continue

        absolute_path = path.locate()
        return str(absolute_path.parent)

    return None


def modify_path():
    """Add the DLL directory to the module search path.

    This will only modify path if
    * on Windows
    * the ``extra-dll`` directory is in package file metadata
    """
    if os.name != "nt":
        return

    # pylint: disable=import-outside-toplevel
    import importlib.metadata

    # pylint: enable=import-outside-toplevel

    try:
        bezier_files = importlib.metadata.files("bezier")
    except importlib.metadata.PackageNotFoundError:
        return

    extra_dll_dir = _get_extra_dll_dir(bezier_files)
    if extra_dll_dir is None:
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
