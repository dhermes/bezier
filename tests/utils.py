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

import os
import pathlib
import shlex
import subprocess
import sys
import textwrap


# See: https://docs.python.org/3/library/platform.html#cross-platform
if sys.maxsize == 2 ** 63 - 1:
    IS_64_BIT = True
elif sys.maxsize == 2 ** 31 - 1:  # pragma: NO COVER
    IS_64_BIT = False
else:  # pragma: NO COVER
    raise ImportError("Unexpected maxsize", sys.maxsize)

IS_MACOS = sys.platform == "darwin"
IS_WINDOWS = os.name == "nt"
IS_LINUX = sys.platform in ("linux", "linux2")
IS_PYPY = sys.implementation.name == "pypy"
INVALID_PATH = "/invalid/path"


def get_gfortran_lib():
    """Find the first directory containing ``libgfortran``.

    This is only intended for Linux or macOS.
    """
    cmd = ("gfortran", "-print-search-dirs")
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    return_code = process.wait()
    if return_code != 0:
        return INVALID_PATH

    cmd_output = process.stdout.read().decode("utf-8")
    parts = cmd_output.split("\nlibraries: =")
    if len(parts) != 2:
        return INVALID_PATH

    library_lines = parts[1].split("\n", 1)
    library_line = library_lines[0]
    # NOTE: ``ctypes.util.find_library()`` can't be used for this because
    #       ``LD_LIBRARY_PATH`` is only set at Python start.
    for part in library_line.split(os.pathsep):
        path = pathlib.Path(part).resolve()
        if list(path.glob("libgfortran*")):
            return str(path)

    return INVALID_PATH


def _invoke_shell(args_str, cmd_directory):
    """Run a command in a shell.

    This is only intended for Linux or macOS.
    """
    args = shlex.split(args_str)
    # NOTE: We print to the stdout of the doctest, rather than using
    #       `subprocess.call()` directly.
    output_bytes = subprocess.check_output(args, cwd=cmd_directory).rstrip()
    print(output_bytes.decode("utf-8"))


def make_invoke_shell(cmd_directory):
    """Produce a function that will run a command in a shell."""

    def wrapped(args_str):
        return _invoke_shell(args_str, cmd_directory)

    return wrapped


def get_git_root():
    """Get the root of the current ``git`` repository."""
    return (
        subprocess.check_output(("git", "rev-parse", "--show-toplevel"))
        .strip()
        .decode("utf-8")
    )


def repo_relative(*path_parts):
    """Get a path relative to the root of the ``git`` repository."""
    git_root = get_git_root()
    return os.path.join(git_root, *path_parts)


def bezier_locate():
    """Locate directories for ``libbezier``.

    In particular, the lib (``-L``) and include (``-I``) directories.
    """
    # NOTE: This will **fail** with a ``KeyError`` if the environment
    #       variable ``BEZIER_INSTALL_PREFIX`` is not set. This is
    #       **intentional**.
    install_prefix = os.environ["BEZIER_INSTALL_PREFIX"]
    # NOTE: This assumes that ``cmake`` (or the build system that installed
    #       ``libbezier``) uses ``include`` and ``lib`` directory names, e.g.
    #       ``/usr/local/include`` and ``/usr/local/lib``.
    bezier_include = os.path.join(install_prefix, "include")
    bezier_lib = os.path.join(install_prefix, "lib")
    if not os.path.isdir(bezier_lib):
        bezier_lib = os.path.join(install_prefix, "lib64")

    return bezier_include, bezier_lib


def _sort_key(name):
    """Sorting helper for members of a directory."""
    return name.lower().lstrip("_")


def tree(directory, suffix=None):
    """Create string (recursively) containing a pretty-printed file tree."""
    names = sorted(os.listdir(directory), key=_sort_key)
    parts = []
    for name in names:
        path = os.path.join(directory, name)
        if os.path.isdir(path):
            sub_part = tree(path, suffix=suffix)
            if sub_part is not None:
                # NOTE: We **always** use posix separator.
                parts.append(name + "/")
                parts.append(textwrap.indent(sub_part, "  "))
        else:
            if suffix is None or name.endswith(suffix):
                if os.path.islink(path):
                    link_dst = os.readlink(path)
                    to_add = f"{name} -> {link_dst}"
                    parts.append(to_add)
                else:
                    parts.append(name)

    if parts:
        return "\n".join(parts)
    else:
        return None


def print_tree(directory, suffix=None):
    """Pretty print a file tree."""
    print(os.path.basename(directory) + os.path.sep)
    full_tree = tree(directory, suffix=suffix)
    print(textwrap.indent(full_tree, "  "))
