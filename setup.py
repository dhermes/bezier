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
import subprocess
import sys

import setuptools
import setuptools.command.build_ext


VERSION = '0.5.0'  # Also in codemeta.json
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
    for part in library_line.split(':'):
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

    This is needed on Mac OS X for a Homebrew installed ``gfortran``, which
    doesn't come with the correct library directories by default (this is
    likely unintentional).

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

    library_dirs = f90_compiler.library_dirs
    # ``library_dirs`` is a list (i.e. mutable), so update in place.
    library_dirs[:] = gfortran_search_path(library_dirs)


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

    patch_library_dirs(f90_compiler)

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
