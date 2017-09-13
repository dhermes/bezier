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

"""Generic (i.e. platform-independent) helpers for ``setup.py``."""

from __future__ import print_function

import collections
import distutils.ccompiler
import os
import subprocess

import setuptools
import setuptools.command.build_ext


GFORTRAN_LIB_ENV = 'GFORTRAN_LIB'
"""Environment variable used to over-ride the ``libgfortran`` search path.

This might be desirable if the "default" ``libgfortran`` does not have some
necessary property. For example, the Homebrew installed ``libgfortran`` on
OS X is not a universal binary. So this environment variable could be used
to point to a custom ``libgfortran`` that **is** universal. (See
``scripts/osx/make_universal_libgfortran.py`` and
``scripts/osx/build-wheels.sh`` for an example of this.)
"""
FORTRAN_LIBRARY_PREFIX = 'libraries: ='
GFORTRAN_MISSING_LIBS = """\
``gfortran`` default library path not found via:

$ gfortran -print-search-dirs
{}"""
GFORTRAN_BAD_PATH = '``gfortran`` library path {} is not a directory.'
BAD_JOURNAL = 'Saving journal failed with {!r}.'
JOURNAL_ENV = 'BEZIER_JOURNAL'
"""Environment variable to specify a text file for saving compiler commands.

Can be used to determine how extension modules were compiled. This can be
useful, for example, to track changes across different systems or simple
to make sure the build is occurring as expected.
"""
# NOTE: This represents the Fortran module dependency graph. Order is
#       important both of the keys and of the dependencies that are in
#       each value.
FORTRAN_MODULES = collections.OrderedDict()
FORTRAN_MODULES['types'] = ('types',)
FORTRAN_MODULES['helpers'] = ('types', 'helpers')
FORTRAN_MODULES['curve'] = ('types', 'helpers', 'curve')
FORTRAN_MODULES['surface'] = ('types', 'helpers', 'curve', 'surface')
FORTRAN_MODULES['curve_intersection'] = (
    'types',
    'helpers',
    'curve',
    'curve_intersection',
)
FORTRAN_SOURCE_FILENAME = os.path.join('src', 'bezier', '{}.f90')
OBJECT_FILENAME = os.path.join('src', 'bezier', '{}.o')
SPEEDUP_FILENAME = os.path.join('src', 'bezier', '_{}_speedup.c')


def gfortran_search_path(library_dirs):
    """Get the library directory paths for ``gfortran``.

    This is a helper for :func:`patch_library_dirs`. Looks for
    ``libraries: =`` in the output of ``gfortran -print-search-dirs`` and
    then parses the paths. If this fails for any reason, this method will
    print an error and return ``library_dirs``.

    Args:
        library_dirs (List[str]): Existing library directories.

    Returns:
        List[str]: The library directories for ``gfortran``.
    """
    cmd = ('gfortran', '-print-search-dirs')
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    return_code = process.wait()

    # Bail out if the command failed.
    if return_code != 0:
        return library_dirs

    cmd_output = process.stdout.read().decode('utf-8')

    # Find single line starting with ``libraries: ``.
    search_lines = cmd_output.strip().split('\n')
    library_lines = [
        line[len(FORTRAN_LIBRARY_PREFIX):]
        for line in search_lines
        if line.startswith(FORTRAN_LIBRARY_PREFIX)
    ]
    if len(library_lines) != 1:
        msg = GFORTRAN_MISSING_LIBS.format(cmd_output)
        print(msg, file=sys.stderr)
        return library_dirs

    # Go through each library in the ``libraries: = ...`` line.
    library_line = library_lines[0]
    accepted = set(library_dirs)
    for part in library_line.split(os.pathsep):
        full_path = os.path.abspath(part.strip())

        if os.path.isdir(full_path):
            accepted.add(full_path)
        else:
            # Ignore anything that isn't a directory.
            msg = GFORTRAN_BAD_PATH.format(full_path)
            print(msg, file=sys.stderr)

    return sorted(accepted)


def patch_library_dirs(f90_compiler):
    """Patch up ``f90_compiler.library_dirs``.

    This assumes (but does not check) that the caller has verified
    ``f90_compiler`` corresponds to ``gfortran``.

    This is needed with ``gfortran`` on OS X and Windows. The
    ``numpy.distutils`` "default" constructor for ``Gnu95FCompiler`` only
    has a single library search path, but there are many library paths
    include in the full ``gcc`` install.

    The ``library_dirs`` directory can be over-ridden by using the
    ``GFORTRAN_LIB`` environment variable (more details on "why" in the
    docstring for (:attr:`GFORTRAN_LIB_ENV`).

    Args:
        f90_compiler (numpy.distutils.fcompiler.gnu.Gnu95FCompiler): A Fortran
            compiler instance.
    """
    gfortran_lib = os.environ.get(GFORTRAN_LIB_ENV)
    # ``library_dirs`` is a list (i.e. mutable), so we can update in place.
    library_dirs = f90_compiler.library_dirs

    if gfortran_lib is None:
        library_dirs[:] = gfortran_search_path(library_dirs)
    else:
        library_dirs[:] = [gfortran_lib]


def extension_modules():
    import numpy as np

    libraries, library_dirs = BuildFortranThenExt.get_library_dirs()
    extensions = []
    for name, dependencies in FORTRAN_MODULES.items():
        if name == 'types':  # No speedup.
            continue

        mod_name = 'bezier._{}_speedup'.format(name)
        path = SPEEDUP_FILENAME.format(name)
        extra_objects = [
            OBJECT_FILENAME.format(dependency)
            for dependency in dependencies
        ]
        extension = setuptools.Extension(
            mod_name,
            [path],
            extra_objects=extra_objects,
            include_dirs=[
                np.get_include(),
                os.path.join('src', 'bezier', 'include'),
            ],
            libraries=libraries,
            library_dirs=library_dirs,
        )
        extensions.append(extension)

    return extensions


class BuildFortranThenExt(setuptools.command.build_ext.build_ext):
    """Custom ``build_ext`` command.

    * Builds Fortran object files for each module (these are required by
      the Cython-generated ``*.c`` extension modules)
    * Provides a static library at ``bezier/lib/libbezier.a`` (on Windows
      this is ``bezier/lib/bezier.lib``)
    * Provides an optional "journaling" feature which allows commands invoked
      during the compilation to be logged (or journaled) to a text file.

    Provides a mutable ``PATCH_FUNCTIONS`` list to allow for "patching" when
    first creating the Fortran compiler ``F90_COMPILER`` that will be
    attached to the class (not the instances).
    """

    # Will be set at runtime, not import time.
    F90_COMPILER = None
    PATCH_FUNCTIONS = []

    def __init__(self, *args, **kwargs):
        setuptools.command.build_ext.build_ext.__init__(self, *args, **kwargs)
        self.journal_file = os.environ.get(JOURNAL_ENV)
        self.commands = []

    @classmethod
    def set_f90_compiler(cls):
        import numpy.distutils.core
        import numpy.distutils.fcompiler

        if cls.F90_COMPILER is not None:
            return

        c_compiler = distutils.ccompiler.new_compiler()
        if c_compiler is None:
            return

        f90_compiler = numpy.distutils.fcompiler.new_fcompiler(
            requiref90=True, c_compiler=c_compiler)
        if f90_compiler is None:
            return

        dist = numpy.distutils.core.get_distribution(always=True)
        f90_compiler.customize(dist)

        # Patch up ``f90_compiler`` with any OS-specific patches.
        for patch_fn in cls.PATCH_FUNCTIONS:
            patch_fn(f90_compiler)

        cls.F90_COMPILER = f90_compiler

    @classmethod
    def has_f90_compiler(cls):
        cls.set_f90_compiler()
        return cls.F90_COMPILER is not None

    @classmethod
    def get_library_dirs(cls):
        cls.set_f90_compiler()
        return cls.F90_COMPILER.libraries, cls.F90_COMPILER.library_dirs

    def start_journaling(self):
        """Capture calls to the system by compilers.

        See: https://github.com/numpy/numpy/blob/v1.13.1/\
        numpy/distutils/ccompiler.py#L154

        Intercepts all calls to ``CCompiler.spawn`` and keeps the
        arguments around to be stored in the local ``commands``
        instance attribute.
        """
        import numpy.distutils.ccompiler

        if self.journal_file is None:
            return

        def journaled_spawn(patched_self, cmd, display=None):
            self.commands.append(cmd)
            return numpy.distutils.ccompiler.CCompiler_spawn(
                patched_self, cmd, display=None)

        numpy.distutils.ccompiler.replace_method(
            distutils.ccompiler.CCompiler,
            'spawn',
            journaled_spawn,
        )

    @staticmethod
    def _command_to_text(command):
        # NOTE: This assumes, but doesn't check that the command has 3
        #       or more arguments.
        first_line = '$ {} \\'
        middle_line = '>   {} \\'
        last_line = '>   {}'

        parts = [first_line.format(command[0])]
        for argument in command[1:-1]:
            parts.append(middle_line.format(argument))
        parts.append(last_line.format(command[-1]))

        return '\n'.join(parts)

    def _commands_to_text(self):
        separator = '-' * 40

        parts = [separator]
        for command in self.commands:
            command_text = self._command_to_text(command)
            parts.extend([command_text, separator])

        parts.append('')  # Trailing newline in file.
        return '\n'.join(parts)

    def save_journal(self):
        """Save journaled commands to file.

        If there is no active journal, does nothing.

        If saving the commands to a file fails, a message will be printed to
        STDERR but the failure will be swallowed so that the extension can
        be built successfully.
        """
        if self.journal_file is None:
            return

        try:
            as_text = self._commands_to_text()
            with open(self.journal_file, 'w') as file_obj:
                file_obj.write(as_text)
        except Exception as exc:
            msg = BAD_JOURNAL.format(exc)
            print(msg, file=sys.stderr)

    def compile_fortran_obj_files(self):
        source_files = [
            FORTRAN_SOURCE_FILENAME.format(mod_name)
            for mod_name in FORTRAN_MODULES
        ]
        obj_files = self.F90_COMPILER.compile(
            source_files,
            output_dir=None,
            macros=[],
            include_dirs=[],
            debug=None,
            extra_postargs=[],
            depends=[],
        )

        # Create a static library (i.e. a ``.a`` / ``.lib`` file).
        c_compiler = self.F90_COMPILER.c_compiler

        static_lib_dir = os.path.join(self.build_lib, 'bezier', 'lib')
        if not os.path.exists(static_lib_dir):
            os.makedirs(static_lib_dir)
        c_compiler.create_static_lib(
            obj_files, 'bezier', output_dir=static_lib_dir)

    def run(self):
        self.set_f90_compiler()
        self.start_journaling()

        self.compile_fortran_obj_files()

        result = setuptools.command.build_ext.build_ext.run(self)
        self.save_journal()
        return result