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
"""Create a universal library for ``libgfortran`` and its dependencies.

Intended to be used on OS X only. It is needed because the Homebrew
installed ``libgfortran`` is ``x86_64`` only, but it also distributes
an ``i386`` version of the necessary dynamic libraries.
"""

from __future__ import print_function

import os
import shutil
import subprocess
import sys

CURR_DIR = os.path.abspath(os.path.dirname(__file__))
FRANKENSTEIN = os.path.join(CURR_DIR, 'frankenstein')
LIBGFORTRAN = 'libgfortran.dylib'
FORTRAN_LIBRARY_PREFIX = 'libraries: ='
LIBRARY_DIRS_ERR = 'Fortran search default library path not found.'


def get_library_dirs():
    """Get library directories in ``gfortran`` search path.

    Uses the information from ``gfortran -print-search-dirs``.

    If the command line output is not of the expected format, prints
    an error message and exits the program with a status code of 1.

    Returns:
        List[str]: Directories in the ``gfortran`` search path.
    """
    cmd = ('gfortran', '-print-search-dirs')
    cmd_output = subprocess.check_output(cmd).decode('utf-8')
    search_lines = cmd_output.strip().split('\n')
    library_lines = [
        line[len(FORTRAN_LIBRARY_PREFIX):]
        for line in search_lines
        if line.startswith(FORTRAN_LIBRARY_PREFIX)
    ]
    if len(library_lines) != 1:
        print(LIBRARY_DIRS_ERR, file=sys.stderr)
        sys.exit(1)
    library_line = library_lines[0]
    directories = []
    for part in library_line.split(':'):
        full_path = os.path.abspath(part)
        if not os.path.exists(full_path):
            continue

        if os.path.isdir(full_path):
            directories.append(full_path)
        else:
            msg = 'Path {} is not a directory.'.format(full_path)
            print(msg, file=sys.stderr)
    if directories:
        print('``gfortran`` library directories:')
        for directory in directories:
            print('\t{}'.format(directory))
    else:
        print('No ``gfortran`` library directories found.', file=sys.stderr)
        sys.exit(1)
    return directories


def find_libgfortran():
    """Get the directory and name of ``libgfortran``.

    Assumes and checks that this ``libgfortran`` is **only** for
    ``x86_64``.

    Exits the program with a status code of 1 if:

    * there is more than one (or zero) directories that contain
      ``libgfortran.dylib``.
    * the ``libgfortran.dylib`` found is **not** ``x86_64``.

    Returns:
        Tuple[str, str]: The directory that contains ``libgfortran`` and the
        full name (with version) of ``libgfortran``.
    """
    library_dirs = get_library_dirs()
    matches = []
    for library_dir in library_dirs:
        path = os.path.join(library_dir, LIBGFORTRAN)
        versioned_path = os.path.realpath(path)
        if os.path.exists(versioned_path):
            matches.append(versioned_path)
    if len(matches) != 1:
        msg = 'Expected exactly one match: {}'.format(', '.join(matches))
        print(msg, file=sys.stderr)
        sys.exit(1)
    dylib = matches[0]
    architectures = get_architectures(dylib)
    if architectures != ['x86_64']:
        msg = 'Expected {} to be x86_64 only, not {}.'.format(
            dylib, ', '.join(architectures)
        )
        print(msg, file=sys.stderr)
        sys.exit(1)
    x86_64_dir, libgfortran = os.path.split(dylib)
    print('Found x86_64 ``libgfortran``:')
    print('\t{}'.format(dylib))
    return x86_64_dir, libgfortran


def get_i386_dir(x86_64_dir, libgfortran):
    """Gets directory containing dynamic libraries targeting ``i386``.

    Exits the program with a status code of 1 if:

    * The expected location of the ``i386`` verison of ``libgfortran``
      does not exist
    * t

    Args:
        x86_64_dir (str): Directory containing ``x86_64`` dynamic libraries.
        libgfortran (str): The name (not path) of the ``libgfortran`` dynamic
            library (should include version).

    Returns:
        str: The directory containing ``i386`` binaries.
    """
    i386_dir = os.path.join(x86_64_dir, 'i386')
    dylib = os.path.join(i386_dir, libgfortran)
    if not os.path.exists(dylib):
        template = 'Expected location of i386 libgfortran does not exist: {}'
        print(template.format(dylib), file=sys.stderr)
        sys.exit(1)
    architectures = get_architectures(dylib)
    if architectures != ['i386']:
        msg = 'Expected {} to be i386 only, not {}.'.format(
            dylib, ', '.join(architectures)
        )
        print(msg, file=sys.stderr)
        sys.exit(1)
    print('Found directory with ``i386`` dynamic libraries:')
    print('\t{}'.format(i386_dir))
    return i386_dir


def get_otool_path(otool_line):
    """Parse path from a line from ``otool -L`` output.

    This **assumes** the format, but does not check it.

    Args:
        otool_line (str): A dependency (or ``install_name``) from ``otool -L``
            output. Expected to be of the form '\t{PATH} (compatibility ...)'.

    Returns:
        str: The ``PATH`` in the ``otool_line``.
    """
    parts = otool_line.split()
    return parts[0]


def get_dependencies(dylib, check_exists=True):
    """Get dependencies for a dynamic library.

    For example:

    .. code-block:: python

       >>> dylib = '/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libgfortran.4.dylib'
       >>> get_dependencies(dylib)
       ['/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libquadmath.0.dylib',
        '/usr/lib/libSystem.B.dylib',
        '/usr/local/lib/gcc/7/libgcc_s.1.dylib']

    Args:
        dylib (str): The path to a dynamic library.
        check_exists (Optional[bool]): Indicates if the existence of a
            dependency should be checked.

    Returns:
        List[str]: The dependencies of ``dylib``.

    Raises:
        ValueError: If the first line of the output is not ``{dylib}:``.
        ValueError: If the ``install_name`` (i.e. the second line) does not
            have the same name (not necessarily same path) as ``dylib``.
            For example
            ``/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libgfortran.4.dylib``
            has an install name of
            ``/usr/local/opt/gcc/lib/gcc/7/libgfortran.4.dylib``.
        ValueError: If one of the dependencies (from any line other than the
            first two) is not actual a file on the current machine. (Will
            only be raised if ``check_exists=True``.)
    """
    cmd = ('otool', '-L', dylib)
    cmd_output = subprocess.check_output(cmd).decode('utf-8')
    lines = cmd_output.strip().split('\n')
    if lines[0] != dylib + ':':
        raise ValueError('Unexpected first line', lines[0])

    install_name = get_otool_path(lines[1])
    if os.path.basename(install_name) != os.path.basename(dylib):
        raise ValueError('Unexpected install_name', install_name, dylib)

    dependencies = []
    for line in lines[2:]:
        dependency = get_otool_path(line)
        if check_exists and not os.path.exists(dependency):
            raise ValueError('Dependency does not exist', dependency)

        dependencies.append(dependency)
    return dependencies


def get_architectures(dylib):
    """Determine the architectures that a dynamic library supports.

    Uses ``lipo -info`` to determine the architectures. Expects outputs to
    resemble one of the following:

    .. code-block:: console

       $ lipo -info /usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libgfortran.dylib
       Non-fat file: .../libgfortran.dylib is architecture: x86_64
       $ lipo -info /usr/local/lib/gcc/7/libgcc_s.1.dylib
       Architectures in the fat file: .../libgcc_s.1.dylib are: x86_64 i386

    (Path information has been replaced by ``...`` for display purposes.)

    Putting this to use:

    .. code-block:: python

       >>> libgfortran = (
       ...     '/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libgfortran.4.dylib')
       >>> get_architectures(libgfortran)
       ['x86_64']
       >>> libgcc_s = '/usr/local/lib/gcc/7/libgcc_s.1.dylib'
       >>> get_architectures(libgcc_s)
       ['x86_64', 'i386']

    Args:
        dylib (str): The path to a dynamic library.

    Returns:
        List[str]: The architecture(s) supported by ``dylib``.

    Raises:
        ValueError: If the ``lipo -info {dylib}`` output does not conform
            to one of the two expected formats.
    """
    cmd = ('lipo', '-info', dylib)
    cmd_output = subprocess.check_output(cmd).decode('utf-8').strip()
    prefix = 'Architectures in the fat file: {} are: '.format(dylib)
    if cmd_output.startswith(prefix):
        architectures = cmd_output[len(prefix):].split()
        return architectures

    else:
        prefix = 'Non-fat file: {} is architecture: '.format(dylib)
        if not cmd_output.startswith(prefix):
            raise ValueError('Unexpected output', cmd_output)

        return [cmd_output[len(prefix):]]


def is_universal(dylib):
    """Checks if a dynamic library is a "universal" binary.

    Uses ``get_architectures`` and checks if both ``i386`` and ``x86_64``
    are supported architectures.

    Args:
        dylib (str): The path to a dynamic library.

    Returns:
        bool: Flag indicating if ``dylib`` is a universal binary.
    """
    architectures = get_architectures(dylib)
    return 'i386' in architectures and 'x86_64' in architectures


def non_universal_libraries(dylib):
    """Get all dependencies (recursively) that are not universal binaries.

    For example:

    .. code-block:: python

       >>> libgfortran = (
       ...     '/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libgfortran.4.dylib')
       >>> non_universal_libraries(libgfortran)
       {'/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libgfortran.4.dylib',
        '/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/libquadmath.0.dylib'}

    Args:
        dylib (str): Path to a dynamic library.

    Returns:
        Set[str]: All non-universal libraries needed by ``dylib``
        (possibly including itself).
    """
    result = set()
    if is_universal(dylib):
        return result

    else:
        result.add(dylib)
    for dependency in get_dependencies(dylib):
        if is_universal(dependency):
            continue

        result.add(dependency)
        result.update(non_universal_libraries(dependency))
    return result


def verify_libraries(libgfortran_path, libraries):
    """Verifies the non-universal dependencies of ``libgfortran``.

    Checks these against our assumption that there are **only** two:
    ``libgfortran`` and ``libquadmath``.

    Exits the program with a status code of 1 if:

    * ``libraries`` does not have two members or ``libgfortran_path`` is
      not one of them.
    * the ``libquadmath`` found is **not** ``x86_64``.
    * the ``libquadmath`` found is not in the same directory as
      ``libgfortran``.

    Args:
        libgfortran_path (str): The full path to the ``libgfortran`` dynamic
            library.
        libraries (Set[str]): The non-universal dynamic libraries required
            by ``libgfortran``. Determined by :func:`non_universal_libraries`.

    Returns:
        str: The name (not path) of the ``libquadmath`` dynamic library.
    """
    print('Non-universal libraries found:')
    for library in libraries:
        print('\t{}'.format(library))
    if len(libraries) != 2 or libgfortran_path not in libraries:
        msg = 'Expected ``libgfortran`` and ``libquadmath``: {}'.format(
            ', '.join(sorted(libraries))
        )
        print(msg, file=sys.stderr)
        sys.exit(1)
    libraries.remove(libgfortran_path)
    libquadmath_path = libraries.pop()
    architectures = get_architectures(libquadmath_path)
    if architectures != ['x86_64']:
        msg = 'Expected {} to be x86_64 only, not {}.'.format(
            libquadmath_path, ', '.join(architectures)
        )
        print(msg, file=sys.stderr)
        sys.exit(1)
    library_dir, libquadmath = os.path.split(libquadmath_path)
    if library_dir != os.path.dirname(libgfortran_path):
        msg = 'Expected {} and {} in same directory.'.format(
            libgfortran_path, libquadmath
        )
        print(msg, file=sys.stderr)
        sys.exit(1)
    return libquadmath


def make_root_dir():
    """Creates "frankenstein" directory for a universal ``libgfortran``.

    If the directory already exists, then exits the program with a status
    code of 1.
    """
    if os.path.exists(FRANKENSTEIN):
        msg = 'The directory {} already exists.'.format(FRANKENSTEIN)
        print(msg, file=sys.stderr)
        sys.exit(1)
    else:
        os.mkdir(FRANKENSTEIN)


def copyfile(source, destination):
    """Copy a file and print a message.

    Args:
        source (str): The file being copied.
        destination (str): The place to copy the file.
    """
    shutil.copyfile(source, destination)
    msg = 'Copied:\n\t  {}\n\t->{}'.format(source, destination)
    print(msg)


def copy_arch(arch, library_dir, libgfortran, libquadmath):
    """Copy libraries specific to a given architecture.

    Args:
        arch (str): The architecture being copied.
        library_dir (str): The directory containing the dynamic libraries.
        libgfortran (str): The name (not path) of the ``libgfortran``
            dynamic library.
        libquadmath (str): The name (not path) of the ``libquadmath``
            dynamic library.

    Returns:
        Tuple[str, str, str, str]: Four-tuple of

        * The path to the ``arch``-specific location of the newly
          created ``libgfortran``
        * The path to the location of the universal ``libgfortran``
          (not yet created, but reference here as the ``install_name``)
        * The path to the ``arch``-specific location of the newly
          created ``libquadmath``
        * The path to the location of the universal ``libquadmath``
          (not yet created, but reference here as the ``install_name``)
    """
    sub_dir = os.path.join(FRANKENSTEIN, arch)
    os.mkdir(sub_dir)
    # Determine the old/new filenames.
    old_libgfortran = os.path.join(library_dir, libgfortran)
    arch_libgfortran = os.path.join(sub_dir, libgfortran)
    universal_libgfortran = os.path.join(FRANKENSTEIN, libgfortran)
    old_libquadmath = os.path.join(library_dir, libquadmath)
    arch_libquadmath = os.path.join(sub_dir, libquadmath)
    universal_libquadmath = os.path.join(FRANKENSTEIN, libquadmath)
    # Update ``libgfortran``
    copyfile(old_libgfortran, arch_libgfortran)
    os.chmod(arch_libgfortran, 0o644)
    subprocess.check_call(
        ('install_name_tool', '-id', universal_libgfortran, arch_libgfortran)
    )
    subprocess.check_call(
        (
            'install_name_tool',
            '-change',
            old_libquadmath,
            universal_libquadmath,
            arch_libgfortran,
        )
    )
    os.chmod(arch_libgfortran, 0o444)
    print('{}:'.format(arch_libgfortran))
    print('\t``install_name``:')
    print('\t\t{}'.format(universal_libgfortran))
    print('\tDependencies:')
    dependencies = get_dependencies(arch_libgfortran, check_exists=False)
    for dependency in dependencies:
        print('\t\t{}'.format(dependency))
    # Update ``libquadmath``
    copyfile(old_libquadmath, arch_libquadmath)
    os.chmod(arch_libquadmath, 0o644)
    subprocess.check_call(
        ('install_name_tool', '-id', universal_libquadmath, arch_libquadmath)
    )
    os.chmod(arch_libquadmath, 0o444)
    print('{}:'.format(arch_libquadmath))
    print('\t``install_name``:')
    print('\t\t{}'.format(universal_libquadmath))
    print('\tDependencies:')
    dependencies = get_dependencies(arch_libquadmath, check_exists=False)
    for dependency in dependencies:
        print('\t\t{}'.format(dependency))
    return (
        arch_libgfortran,
        universal_libgfortran,
        arch_libquadmath,
        universal_libquadmath,
    )


def combine_dylibs(i386_dylib, x86_64_dylib, universal_dylib):
    """Combine two dynamic libraries into one universal dynamic library.

    Args:
        i386_dylib (str): The full path to the copy of the dynamic library
            that targets the ``i386`` architecture.
        x86_64_dylib (str): The full path to the copy of the dynamic library
            that targets the ``x86_64`` architecture.
        universal_dylib (str): The full path of the universal dynamic library
            that they should be combined into.
    """
    subprocess.check_call(
        (
            'lipo',
            i386_dylib,
            x86_64_dylib,
            '-create',
            '-output',
            universal_dylib,
        )
    )
    print('Created universal dynamic library:')
    print('\t{}'.format(universal_dylib))
    curr_dir = os.getcwd()
    # Make a symlink **without** the library version.
    # NOTE: This assumes that os.path.dirname(universal_dylib) == FRANKENSTEIN.
    filename = os.path.basename(universal_dylib)
    name, _, extension = filename.split('.')
    unversioned = '{}.{}'.format(name, extension)
    os.chdir(FRANKENSTEIN)
    os.symlink(filename, unversioned)
    os.chdir(curr_dir)
    print('Created symbolic link:')
    print('\t{}@ -> {}'.format(unversioned, filename))


def main():
    make_root_dir()
    x86_64_dir, libgfortran = find_libgfortran()
    i386_dir = get_i386_dir(x86_64_dir, libgfortran)
    full_path = os.path.join(x86_64_dir, libgfortran)
    libraries = non_universal_libraries(full_path)
    libquadmath = verify_libraries(full_path, libraries)
    i386_names = copy_arch('i386', i386_dir, libgfortran, libquadmath)
    x86_64_names = copy_arch('x86_64', x86_64_dir, libgfortran, libquadmath)
    (
        i386_libgfortran,
        universal_libgfortran,
        i386_libquadmath,
        universal_libquadmath,
    ) = i386_names
    x86_64_libgfortran, _, x86_64_libquadmath, _ = x86_64_names
    combine_dylibs(i386_libgfortran, x86_64_libgfortran, universal_libgfortran)
    combine_dylibs(i386_libquadmath, x86_64_libquadmath, universal_libquadmath)


if __name__ == '__main__':
    main()
