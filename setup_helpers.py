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
import copy
import distutils.ccompiler
import os
import shutil
import subprocess
import sys

import setuptools
import setuptools.command.build_ext


DEBUG_ENV = 'DEBUG'
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
# NOTE: These are mostly recommendations from Certik found here:
#         http://www.fortran90.org/src/faq.html
#       specifically "What compiler options should I use for ...?".
FPIC = '-fPIC'
GFORTRAN_SHARED_FLAGS = (  # Used for both "DEBUG" and "OPTIMIZE"
    '-Wall',
    '-Wextra',
    # ``-Wextra`` includes ``no-compare-reals``, which warns about
    # ``value == 0.0_dp``
    '-Wno-compare-reals',
    '-Wimplicit-interface',
    FPIC,
    '-fmax-errors=1',
    '-std=f2008',
)
GFORTRAN_DEBUG_FLAGS = (
    '-g',
    '-fcheck=all',
    '-fbacktrace',
    '-fimplicit-none',
    '-pedantic',
)
GFORTRAN_OPTIMIZE_FLAGS = (
    '-Werror',
    '-O3',
    '-funroll-loops',
)
BAD_JOURNAL = 'Saving journal failed with {!r}.'
JOURNAL_ENV = 'BEZIER_JOURNAL'
"""Environment variable to specify a text file for saving compiler commands.

Can be used to determine how extension modules were compiled. This can be
useful, for example, to track changes across different systems or simple
to make sure the build is occurring as expected.
"""
QUADPACK_DIR = 'quadpack'
QUADPACK_SOURCE_FILENAME = os.path.join('src', 'bezier', QUADPACK_DIR, '{}.f')
# NOTE: QUADPACK module dependencies: order is important.
QUADPACK_MODULES = (
    'd1mach',
    'dqelg',
    'dqpsrt',
    'dqk21',
    'dqagse',
)
# NOTE: This represents the Fortran module dependency graph. Order is
#       important both of the keys and of the dependencies that are in
#       each value.
FORTRAN_MODULES = (
    'types',
    'status',
    'helpers',
    'curve',
    'surface',
    'curve_intersection',
    'surface_intersection',
)
FORTRAN_SOURCE_FILENAME = os.path.join('src', 'bezier', '{}.f90')
OBJECT_FILENAME = os.path.join('src', 'bezier', '{}.o')
SPEEDUP_FILENAME = os.path.join('src', 'bezier', '_speedup.c')


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
    if BuildFortranThenExt.USE_SHARED_LIBRARY:
        # Here we don't depend on object files since the functionality
        # is contained in the shared library.
        extra_objects = []
    else:
        # NOTE: These may be treated as relative paths and replaced
        #       before the extension is actually built.
        extra_objects = [
            OBJECT_FILENAME.format(fortran_module)
            for fortran_module in FORTRAN_MODULES
        ]
        extra_objects.extend(
            OBJECT_FILENAME.format(os.path.join(QUADPACK_DIR, fortran_module))
            for fortran_module in QUADPACK_MODULES
        )

    # NOTE: Copy ``libraries`` and ``library_dirs`` so they
    #       aren't shared and mutable.
    extension = setuptools.Extension(
        'bezier._speedup',
        [SPEEDUP_FILENAME],
        extra_objects=extra_objects,
        include_dirs=[
            np.get_include(),
            os.path.join('src', 'bezier', 'include'),
        ],
        libraries=copy.deepcopy(libraries),
        library_dirs=copy.deepcopy(library_dirs),
    )

    return [extension]


def patch_f90_compiler(f90_compiler):
    """Patch up ``f90_compiler``.

    For now, only updates the flags for ``gfortran``. In this case, it add
    any of ``GFORTRAN_SHARED_FLAGS`` that are missing. In debug mode, it also
    adds any flags in ``GFORTRAN_DEBUG_FLAGS`` and makes sure none of the flags
    in ``GFORTRAN_OPTIMIZE_FLAGS`` are present. In standard mode ("OPTIMIZE"),
    makes sure flags in ``GFORTRAN_OPTIMIZE_FLAGS`` are present and flags in
    ``GFORTRAN_DEBUG_FLAGS`` are not.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.
    """
    # NOTE: NumPy may not be installed, but we don't want **this** module to
    #       cause an import failure in ``setup.py``.
    from numpy.distutils.fcompiler import gnu

    # Only ``gfortran``.
    if not isinstance(f90_compiler, gnu.Gnu95FCompiler):
        return False

    # NOTE: Should probably handle ``compiler_f77`` too.
    f90_flags = f90_compiler.compiler_f90
    for flag in GFORTRAN_SHARED_FLAGS:
        if flag not in f90_flags:
            f90_flags.append(flag)

    if DEBUG_ENV in os.environ:
        to_add = GFORTRAN_DEBUG_FLAGS
        to_remove = GFORTRAN_OPTIMIZE_FLAGS
    else:
        to_add = GFORTRAN_OPTIMIZE_FLAGS
        to_remove = GFORTRAN_DEBUG_FLAGS

    for flag in to_add:
        if flag not in f90_flags:
            f90_flags.append(flag)

    without = [flag for flag in f90_flags
               if flag not in to_remove]
    # Update in place.
    f90_flags[:] = without


class BuildFortranThenExt(setuptools.command.build_ext.build_ext):
    """Custom ``build_ext`` command.

    * Builds Fortran object files for each module (these are required by
      the Cython-generated ``*.c`` extension modules)
    * Provides a static library at ``bezier/lib/libbezier.a`` (on Windows
      this is ``bezier/lib/bezier.lib``)
    * Provides an optional "journaling" feature which allows commands invoked
      during the compilation to be logged (or journaled) to a text file.

    Provides mutable class attributes:

    * ``PATCH_FUNCTIONS`` list to allow for "patching" when first creating
      the Fortran compiler ``F90_COMPILER`` that will be attached to the
      class (not the instances).
    * ``CUSTOM_STATIC_LIB`` callable that takes a list of object files and
      uses them to create a static / shared library. If not provided, then
      :meth:`_default_static_lib` will be used.
    * ``USE_SHARED_LIBRARY`` flag indicating if extensions will contain built
      object files or if they will refer to a shared library.
    * ``CLEANUP`` optional callable that cleans up at the end of :meth:`run`.
    """

    # Will be set at runtime, not import time.
    F90_COMPILER = None
    PATCH_FUNCTIONS = []
    CUSTOM_STATIC_LIB = None
    USE_SHARED_LIBRARY = False
    CLEANUP = None

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
        if c_compiler.compiler_type == 'msvc':
            c_compiler.initialize()

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
        if cls.USE_SHARED_LIBRARY:
            # NOTE: This assumes that the `libbezier` shared library will
            #       contain all libraries needed (e.g. there is no
            #       dependendence on ``libgfortran`` or similar). It's expected
            #       that ``library_dirs`` will be updated at run-time to have
            #       temporary build directories added.
            libraries = ['bezier']
            library_dirs = []
            return libraries, library_dirs
        else:
            return cls.F90_COMPILER.libraries, cls.F90_COMPILER.library_dirs

    def start_journaling(self):
        """Capture calls to the system by compilers.

        See: https://github.com/numpy/numpy/blob/v1.14.0/\
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

    def _default_static_lib(self, obj_files):
        """Create a static library (i.e. a ``.a`` / ``.lib`` file).

        Args:
            obj_files (List[str]): List of paths of compiled object files.
        """
        c_compiler = self.F90_COMPILER.c_compiler

        static_lib_dir = os.path.join(self.build_lib, 'bezier', 'lib')
        if not os.path.exists(static_lib_dir):
            os.makedirs(static_lib_dir)
        c_compiler.create_static_lib(
            obj_files, 'bezier', output_dir=static_lib_dir)

        # NOTE: We must "modify" the paths for the ``extra_objects`` in
        #       each extension since they were compiled with
        #       ``output_dir=self.build_temp``.
        for extension in self.extensions:
            extension.extra_objects[:] = [
                os.path.join(self.build_temp, rel_path)
                for rel_path in extension.extra_objects
            ]

    def compile_fortran_obj_files(self):
        source_files_quadpack = [
            QUADPACK_SOURCE_FILENAME.format(mod_name)
            for mod_name in QUADPACK_MODULES
        ]
        source_files_bezier = [
            FORTRAN_SOURCE_FILENAME.format(mod_name)
            for mod_name in FORTRAN_MODULES
        ]
        source_files = source_files_quadpack + source_files_bezier
        obj_files = self.F90_COMPILER.compile(
            source_files,
            output_dir=self.build_temp,
            macros=[],
            include_dirs=[],
            debug=None,
            extra_postargs=[
                '-J',
                self.build_temp,
            ],
            depends=[],
        )

        if self.CUSTOM_STATIC_LIB is None:
            self._default_static_lib(obj_files)
        else:
            self.CUSTOM_STATIC_LIB(obj_files)

    def _default_cleanup(self):
        """Default cleanup after :meth:`run`.

        For in-place builds, moves the built shared library into the source
        directory.
        """
        if not self.inplace:
            return

        shutil.move(
            os.path.join(self.build_lib, 'bezier', 'lib'),
            os.path.join('src', 'bezier'),
        )

    def cleanup(self):
        """Cleanup after :meth:`run`."""
        if self.CLEANUP is None:
            self._default_cleanup()
        else:
            self.CLEANUP()

    def run(self):
        self.set_f90_compiler()
        self.start_journaling()

        self.compile_fortran_obj_files()

        result = setuptools.command.build_ext.build_ext.run(self)
        self.save_journal()
        self.cleanup()

        return result
