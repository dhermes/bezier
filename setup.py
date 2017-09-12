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

"""Setup file for bezier."""

from __future__ import print_function

import collections
import distutils.ccompiler
import os
import pkg_resources
import platform
import shutil
import subprocess
import sys
import tempfile

import setuptools
import setuptools.command.build_ext


VERSION = '0.5.0.dev1'  # Also in codemeta.json
PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))
README_FILENAME = os.path.join(PACKAGE_ROOT, 'README.rst')
NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ pip              install numpy
$ python    -m pip install numpy --user
$ python2.7 -m pip install numpy --user
$ python3.6 -m pip install numpy --user
$ # OR
$ conda install numpy
"""
MISSING_F90_MESSAGE = """\
No Fortran 90 compiler found.

Skipping Fortran extension speedups.
"""
WINDOWS_MESSAGE = """\
Skipping Fortran extension speedups on Windows.

Sorry for this inconvenience. For more information or to help, visit:

    https://github.com/dhermes/bezier/issues/26
"""
REQUIREMENTS = (
    'numpy >= 1.11.2',
    'six >= 1.9.0',
)
EXTRAS_REQUIRE = {
    ':python_version<"3.4"': ['enum34'],
}
DESCRIPTION = (
    u'Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects')
FORTRAN_LIBRARY_PREFIX = 'libraries: ='
GFORTRAN_ERR_MSG = '``gfortran`` default library path not found.'
GFORTRAN_BAD_PATH = '``gfortran`` library path {} is not a directory.'
BAD_JOURNAL = 'Saving journal failed with {!r}.'
JOURNAL_ENV = 'BEZIER_JOURNAL'
GFORTRAN_LIB_ENV = 'GFORTRAN_LIB'
MAC_OS_X = 'darwin'
# NOTE: This represents the Fortran module dependency graph. Order is
#       important both of the keys and of the dependencies that are in
#       each value.
FORTRAN_MODULES = collections.OrderedDict()
FORTRAN_MODULES['types'] = ('types',)
FORTRAN_MODULES['helpers'] = ('types', 'helpers')
FORTRAN_MODULES['curve'] = ('types', 'curve')
FORTRAN_MODULES['surface'] = ('types', 'curve', 'surface')
FORTRAN_MODULES['curve_intersection'] = (
    'types',
    'helpers',
    'curve',
    'curve_intersection',
)
FORTRAN_SOURCE_FILENAME = os.path.join('src', 'bezier', '{}.f90')
OBJECT_FILENAME = os.path.join('src', 'bezier', '{}.o')
SPEEDUP_FILENAME = os.path.join('src', 'bezier', '_{}_speedup.c')


def is_installed(requirement):
    try:
        pkg_resources.require(requirement)
    except pkg_resources.ResolutionError:
        return False
    else:
        return True


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
        print(GFORTRAN_ERR_MSG, file=sys.stderr)
        return library_dirs

    # Go through each library in the ``libraries: = ...`` line.
    library_line = library_lines[0]
    accepted = set(library_dirs)
    for part in library_line.split(os.pathsep):
        full_path = os.path.abspath(part.strip())

        if os.path.isdir(full_path):
            accepted.add(full_path)
        else:
            # Ignore anything that isn't a path.
            msg = GFORTRAN_BAD_PATH.format(full_path)
            print(msg, file=sys.stderr)

    return sorted(accepted)


def patch_library_dirs(f90_compiler):
    """Patch up ``f90_compiler.library_dirs``.

    Assumes the caller is :func:`patch_osx_gfortran`.

    This is needed on Mac OS X for a Homebrew installed ``gfortran``, which
    doesn't come with the correct library directories by default (this is
    likely unintentional).

    The ``library_dirs`` directory can be over-ridden by using the
    ``GFORTRAN_LIB`` environment variable. This might be desirable, since the
    Homebrew ``libgfortran`` is **also** not a universal binary. So this
    over-ride could be used to point to a custom made ``libgfortran.dylib``
    that is a combined version of the ``i386`` and ``x86_64`` versions of
    ``libgfortran`` provided by Homebrew.

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


def is_dual_architecture():
    """Checks if the current Python binary is dual architecture.

    Expected to be used by :func:`patch__compile`.

    Only relevant on OS X. This uses ``lipo -info`` to check that the
    executable is a "fat file" with both ``i386`` and ``x86_64``
    architectures.

    We use ``lipo -info`` rather than ``file`` because ``lipo`` is
    purpose-built for checking the architecture(s) in a file.

    This property could also be checked by looking for the presence of
    multiple architectures in
    ``distutils.sysconfig.get_config_var('LDFLAGS')``.

    Returns:
        bool: Indicating if the Python binary is dual architecture
        (:data:`True`) or single architecture (:data:`False`).
    """
    if sys.platform != MAC_OS_X:
        return False

    cmd = ('lipo', '-info', sys.executable)
    cmd_output = subprocess.check_output(cmd).decode('utf-8').strip()

    prefix = 'Architectures in the fat file: {} are: '.format(sys.executable)

    if cmd_output.startswith(prefix):
        architectures = cmd_output[len(prefix):].split()
        return 'i386' in architectures and 'x86_64' in architectures
    else:
        return False


def gfortran_supports_dual_architecture():
    """Simple check if ``gfortran`` supports dual architecture.

    Only relevant on OS X. By default, the Homebrew ``gfortran`` **does not**
    support building dual architecture object files. This checks support
    for this feature by trying to build a very simple Fortran 90 program with
    ``-arch i386 -arch x86_64``.

    Returns:
        bool: Indicates if ``gfortran`` can build dual architecture binaries.
    """
    if sys.platform != MAC_OS_X:
        return False

    temporary_directory = tempfile.mkdtemp(suffix='-fortran')
    source_name = os.path.join(temporary_directory, 'bar.f90')
    with open(source_name, 'w') as file_obj:
        file_obj.writelines([
            'subroutine bar(x, y)\n',
            '  integer, intent(in) :: x\n',
            '  integer, intent(out) :: y\n',
            '\n',
            '  y = x + 2\n',
            '\n',
            'end subroutine bar\n',
            '\n',
        ])

    object_name = os.path.join(temporary_directory, 'bar.o')
    cmd = (
        'gfortran',
        '-arch', 'i386', '-arch', 'x86_64',
        '-c', source_name,
        '-o', object_name,
    )

    cmd_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    result = b'arch flags ignored' not in cmd_output

    # Clean-up the temporary directory (could be done with a context manager).
    shutil.rmtree(temporary_directory)

    return result


class _DualArchitectureCompile(object):
    """Callable wrapper that over-rides ``_compile``.

    Intended to be used by :func:`patch__compile`.

    Only relevant on OS X. Objects of this type are intended to be used
    to replace / augment to ``_compile`` method on a ``Gnu95FCompiler`` (i.e.
    a Fortran compiler that uses ``gfortran``). This is because the Homebrew
    ``gfortran`` can't build fat binaries:

    .. code-block:: console

       $ gfortran -arch i386 -arch x86_64 -c bar.f90 -o bar.o
       gfortran: warning: x86_64 conflicts with i386 (arch flags ignored)

    So instead, this will compile two separate object files and combine them:

    .. code-block:: console

       $ gfortran -arch i386   -c bar.f90 -o ${I386_DIR}/bar.o
       $ gfortran -arch x86_64 -c bar.f90 -o ${X86_64_DIR}/bar.o
       $ lipo ${I386_DIR}/bar.o ${X86_64_DIR}/bar.o -create -output bar.o

    Args:
        f90_compiler (numpy.distutils.fcompiler.gnu.Gnu95FCompiler): A Fortran
            compiler instance corresponding to ``gfortran``.
    """

    def __init__(self, f90_compiler):
        self.f90_compiler = f90_compiler
        self.original_compile = f90_compiler._compile
        self.compiler_cmd = f90_compiler.compiler_f90
        self.arch_index = None  # Set in ``_verify()``.
        self.arch_value = None  # Set in ``_verify()``.
        self._verify()

    def _verify(self):
        """Makes sure the constructor arguments are valid.

        In particular, makes sure that ``compiler_cmd`` has exactly one instance
        of ``-arch``.

        If this succeeds, will set ``arch_index`` and ``arch_value`` on
        the instance.

        Raises:
            TypeError: If ``compiler_cmd`` is not a ``list``.
            ValueError: If ``compiler_cmd`` doesn't have exactly one ``-arch``
                segment.
            ValueError: If ``-arch`` is the **last** segment in
                ``compiler_cmd``.
            ValueError: If the ``-arch`` value is not ``i386`` or ``x86_64``.
        """
        if not isinstance(self.compiler_cmd, list):
            raise TypeError('Expected a list', self.compiler_cmd)

        if self.compiler_cmd.count('-arch') != 1:
            raise ValueError(
                'Did not find exactly one "-arch" in', self.compiler_cmd)

        arch_index = self.compiler_cmd.index('-arch') + 1
        if arch_index == len(self.compiler_cmd):
            raise ValueError(
                'There is no architecture specified in', self.compiler_cmd)

        arch_value = self.compiler_cmd[arch_index]
        if arch_value not in ('i386', 'x86_64'):
            raise ValueError(
                'Unexpected architecture', arch_value, 'in', self.compiler_cmd)

        self.arch_index = arch_index
        self.arch_value = arch_value

    def _set_architecture(self, architecture):
        """Set the architecture on the Fortran compiler.

        ``compiler_cmd`` is actually a list (mutable), so we can update it here
        and it will change the architecture that ``f90_compiler`` targets.

        Args:
            architecture (str): One of ``i386`` or ``x86_64``.
        """
        self.compiler_cmd[self.arch_index] = architecture

    def _restore_architecture(self):
        """Restore the architecture on the Fortran compiler.

        Resets the ``-arch`` value in ``compiler_cmd`` to its original value.
        """
        self.compiler_cmd[self.arch_index] = self.arch_value

    def __call__(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
        """Call-able replacement for ``_compile``.

        This assumes (but does not verify) that ``original_compile`` has
        no return value.

        Args:
            obj (str): The location of the object file to be created.
            src (str): The location of the source file to be compiled.
            ext (str): The file extension (used to determine flags).
            cc_args (List[str]): Compile args, typically just ``['-c']``.
            extra_postargs (List[str]): Extra arguments at the end of the
                compile command.
            pp_opts (List[str]): Unused by the NumPy ``distutils`` Fortran
                compilers. List of pre-processor options.
        """
        obj_name = os.path.basename(obj)

        # Create a directory and compile an object targeting i386.
        i386_dir = tempfile.mkdtemp(suffix='-i386')
        i386_obj = os.path.join(i386_dir, obj_name)
        self._set_architecture('i386')
        self.original_compile(
            i386_obj, src, ext, cc_args, extra_postargs, pp_opts)

        # Create a directory and compile an object targeting x86_64.
        x86_64_dir = tempfile.mkdtemp(suffix='-x86_64')
        x86_64_obj = os.path.join(x86_64_dir, obj_name)
        self._set_architecture('x86_64')
        self.original_compile(
            x86_64_obj, src, ext, cc_args, extra_postargs, pp_opts)

        # Restore the compiler back to how it was before we modified it (could
        # be done with a context manager).
        self._restore_architecture()

        # Use ``lipo`` to combine the object files into a universal.
        lipo_cmd = ('lipo', i386_obj, x86_64_obj, '-create', '-output', obj)
        self.f90_compiler.spawn(lipo_cmd)

        # Clean up the temporary directories (could be done with a context
        # manager).
        shutil.rmtree(i386_dir)
        shutil.rmtree(x86_64_dir)


def patch__compile(f90_compiler):
    """Modify the Fortran compiler to create universal binary object files.

    Assumes the caller is :func:`patch_osx_gfortran`.

    Does so by patching ``f90_compiler._compile`` with a custom command.

    Patching is only done if the current Python is a universal binary (i.e.
    dual architecture) and the installed ``gfortran`` cannot create universal
    (sometimes called "fat") binaries.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.
    """
    # Only if Python is a universal binary.
    if not is_dual_architecture():
        return

    # Only if ``gfortran`` can't produce universal binaries.
    if gfortran_supports_dual_architecture():
        return

    f90_compiler._compile = _DualArchitectureCompile(f90_compiler)


def patch_osx_gfortran(f90_compiler):
    """Patch up ``f90_compiler.library_dirs``.

    On Mac OS X for a Homebrew installed ``gfortran`` has a few quirks that
    must be dealt with.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.
    """
    from numpy.distutils.fcompiler import gnu

    # Only Mac OS X.
    if sys.platform != MAC_OS_X:
        return
    # Only ``gfortran``.
    if not isinstance(f90_compiler, gnu.Gnu95FCompiler):
        return

    patch_library_dirs(f90_compiler)
    patch__compile(f90_compiler)


def get_f90_compiler():
    import numpy.distutils.core
    import numpy.distutils.fcompiler

    c_compiler = distutils.ccompiler.new_compiler()
    if c_compiler is None:
        return None

    f90_compiler = numpy.distutils.fcompiler.new_fcompiler(
        requiref90=True, c_compiler=c_compiler)
    if f90_compiler is None:
        return None

    dist = numpy.distutils.core.get_distribution(always=True)
    f90_compiler.customize(dist)

    patch_osx_gfortran(f90_compiler)

    return f90_compiler


def compile_fortran_obj_files(f90_compiler, static_lib_dir):
    source_files = [
        FORTRAN_SOURCE_FILENAME.format(mod_name)
        for mod_name in FORTRAN_MODULES
    ]
    obj_files = f90_compiler.compile(
        source_files,
        output_dir=None,
        macros=[],
        include_dirs=[],
        debug=None,
        extra_postargs=[],
        depends=[],
    )

    # Create a .a/.lib file (i.e. a static library).
    c_compiler = f90_compiler.c_compiler
    if not os.path.exists(static_lib_dir):
        os.makedirs(static_lib_dir)
    c_compiler.create_static_lib(
        obj_files, 'bezier', output_dir=static_lib_dir)


def _extension_modules():
    # NOTE: This assumes is_installed('numpy') has already passed.
    #       H/T to https://stackoverflow.com/a/41575848/1068170
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


def extension_modules():
    if platform.system().lower() == 'windows':
        print(WINDOWS_MESSAGE, file=sys.stderr)
        return []
    elif BuildFortranThenExt.has_f90_compiler():
        return _extension_modules()
    else:
        print(MISSING_F90_MESSAGE, file=sys.stderr)
        return []


def make_readme():
    with open(README_FILENAME, 'r') as file_obj:
        return file_obj.read()


class BuildFortranThenExt(setuptools.command.build_ext.build_ext):

    # Will be set at runtime, not import time.
    F90_COMPILER = None

    def __init__(self, *args, **kwargs):
        setuptools.command.build_ext.build_ext.__init__(self, *args, **kwargs)
        self.journal_file = os.environ.get(JOURNAL_ENV)
        self.commands = []

    @classmethod
    def set_f90_compiler(cls):
        if cls.F90_COMPILER is None:
            cls.F90_COMPILER = get_f90_compiler()

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

    def run(self):
        self.set_f90_compiler()
        self.start_journaling()

        static_lib_dir = os.path.join(self.build_lib, 'bezier', 'lib')
        compile_fortran_obj_files(self.F90_COMPILER, static_lib_dir)

        result = setuptools.command.build_ext.build_ext.run(self)
        self.save_journal()
        return result


def setup():
    setuptools.setup(
        name='bezier',
        version=VERSION,
        description=DESCRIPTION,
        author='Danny Hermes',
        author_email='daniel.j.hermes@gmail.com',
        long_description=make_readme(),
        scripts=(),
        url='https://github.com/dhermes/bezier',
        packages=['bezier'],
        package_dir={'': 'src'},
        license='Apache 2.0',
        platforms='Posix; MacOS X; Windows',
        package_data={
            'bezier': [
                '*.pxd',
                os.path.join('include', '*.h'),
                os.path.join('include', 'bezier', '*.h'),
                os.path.join('lib', '*.a'),
                os.path.join('lib', '*.lib'),
            ],
        },
        zip_safe=True,
        install_requires=REQUIREMENTS,
        extras_require=EXTRAS_REQUIRE,
        ext_modules=extension_modules(),
        classifiers=(
            'Development Status :: 4 - Beta',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Mathematics',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ),
        cmdclass={'build_ext': BuildFortranThenExt},
    )


def main():
    if not is_installed('numpy>=1.9.0'):
        print(NUMPY_MESSAGE, file=sys.stderr)
        sys.exit(1)

    setup()


if __name__ == '__main__':
    main()
