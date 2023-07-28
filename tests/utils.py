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

import functools
import os
import pathlib
import re
import shlex
import subprocess
import sys
import textwrap


# See: https://docs.python.org/3/library/platform.html#cross-platform
if sys.maxsize == 2**63 - 1:
    IS_64_BIT = True
elif sys.maxsize == 2**31 - 1:  # pragma: NO COVER
    IS_64_BIT = False
else:  # pragma: NO COVER
    raise ImportError("Unexpected maxsize", sys.maxsize)

IS_MACOS = sys.platform == "darwin"
IS_WINDOWS = os.name == "nt"
IS_LINUX = sys.platform in ("linux", "linux2")
IS_PYPY = sys.implementation.name == "pypy"


def _invoke_shell(args_str, cmd_directory):
    """Run a command in a shell.

    This is only intended for Linux or macOS.
    """
    args = shlex.split(args_str)
    # NOTE: We print to the stdout of the doctest, rather than using
    #       `subprocess.call()` directly.
    output_bytes = subprocess.check_output(args, cwd=cmd_directory)
    print(output_bytes.decode("utf-8"), end="")


def make_invoke_shell(cmd_directory):
    """Produce a function that will run a command in a shell."""

    def wrapped(args_str):
        return _invoke_shell(args_str, cmd_directory)

    return wrapped


def _strip_shell(cmd):
    """Strip a (potentially multi-line) command.

    For example

    $ foo \
    >     --bar baz \
    >     --quux 10

    would becomes `foo --bar baz --quux 10`.
    """
    parts = []
    for line in cmd.strip().split("\n"):
        without_console = line.lstrip("$> ").rstrip("\\ ")
        parts.append(without_console)
    return " ".join(parts)


def _gcc_homebrew_version(gcc_root, gcc_bin):
    """Determine a semver version for a Homebrew installed ``gcc``.

    This includes **four** parts: major, minor, patch and revision. For
    example, ``/usr/local/Cellar/gcc/10.2.0/bin/gcc-10`` or
    ``/usr/local/Cellar/gcc/10.2.0_3/bin/gcc-10`` are examples without and with
    a ``_N`` revision.

    This function does no error handling; the expectation is that the caller
    wraps any error with helpful information e.g. the arguments.
    """
    relative_path = gcc_bin.relative_to(gcc_root)
    version_dir = relative_path.parts[0]

    revision = 0
    if "_" in version_dir:
        # This also asserts length 2.
        version_dir, revision_str = version_dir.split("_")
        revision = int(revision_str)

    major_str, minor_str, patch_str = version_dir.split(".")
    return int(major_str), int(minor_str), int(patch_str), revision


def _find_gcc_homebrew():
    gcc_root_bytes = subprocess.check_output(("brew", "--cellar", "gcc"))
    gcc_root = gcc_root_bytes.decode("utf-8").rstrip()
    matches = list(pathlib.Path(gcc_root).glob("*/bin/gcc-[0-9]*"))

    if len(matches) == 0:
        raise ValueError(
            "Could not find Homebrew-installed ``gcc``",
            gcc_root,
        )

    if len(matches) == 1:
        return str(matches[0])

    sort_func = functools.partial(_gcc_homebrew_version, gcc_root)
    matches.sort(key=sort_func)
    chosen = str(matches[-1])

    matches_str = ", ".join(str(match) for match in matches)
    print(
        f"Found multiple matches ({matches_str}) for Homebrew-installed "
        f"``gcc``, using the newest one: {chosen}.",
        file=sys.stderr,
    )
    return chosen


def _find_gcc():
    if IS_LINUX:
        return "gcc"

    if IS_MACOS:
        return _find_gcc_homebrew()

    if IS_WINDOWS:
        raise NotImplementedError

    raise OSError("Unexpected operating system")


def build_and_run_c(filename):
    """Build and run a C example from ``docs/abi/``."""
    bezier_include, bezier_lib = bezier_locate()
    print(f"$ INCLUDE_DIR={bezier_include}")
    print(f"$ LIB_DIR={bezier_lib}")

    docs_abi_directory = repo_relative("docs", "abi")
    invoke_shell = make_invoke_shell(docs_abi_directory)

    build_pretty = "\n".join(
        [
            "$ gcc \\",
            ">     -o example \\",
            f">     {filename} \\",
            '>     -I "${INCLUDE_DIR}" \\',
            '>     -L "${LIB_DIR}" \\',
            '>     -Wl,-rpath,"${LIB_DIR}" \\',
            ">     -lbezier \\",
            ">     -lm -lgfortran",
        ]
    )
    print(build_pretty)
    gcc_bin = _find_gcc()
    build_pretty = (
        build_pretty.replace("$ gcc", f"$ {gcc_bin}")
        .replace("${INCLUDE_DIR}", bezier_include)
        .replace("${LIB_DIR}", bezier_lib)
    )
    invoke_shell(_strip_shell(build_pretty))

    run_pretty = "$ ./example"
    print(run_pretty)
    invoke_shell(_strip_shell(run_pretty))

    os.unlink(os.path.join(docs_abi_directory, "example"))


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
                parts.append(name + os.path.sep)
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


def print_tree(directory, suffix=None, replacements=None):
    """Pretty print a file tree."""
    if replacements is None:
        replacements = ()

    full_tree = tree(directory, suffix=suffix)
    content = "\n".join(
        [
            os.path.basename(directory) + os.path.sep,
            textwrap.indent(full_tree, "  "),
        ]
    )

    for pattern, replacement in replacements:
        output_str = re.sub(pattern, replacement, output_str)

    print(content)
