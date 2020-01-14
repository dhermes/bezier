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

"""Helpers for running ``doctest``.

This is to enable re-use and to avoid having large amounts of code in
``testsetup`` blocks. This will explicitly not be unit tested and so code
elements may stay around even after there are no more users (without more
diligent checks).
"""

import os
import pathlib
import shlex
import subprocess


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
