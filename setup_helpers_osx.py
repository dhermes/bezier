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

"""Helpers for ``setup.py`` specific to OS X."""

import os
import shutil
import subprocess
import sys
import tempfile

import setup_helpers


MAC_OS_X = 'darwin'
# NOTE: Used by ``gfortran_supports_dual_architecture()``.
SIMPLE_F90_SOURCE = """\
subroutine bar(x, y)
  integer, intent(in) :: x
  integer, intent(out) :: y

  y = x + 2

end subroutine bar
"""


def is_osx_gfortran(f90_compiler):
    """Checks if the current build is ``gfortran`` on OS X.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.

    Returns:
        bool: Only :data:`True` if

        * Current OS is OS X (checked via ``sys.platform``).
        * ``f90_compiler`` corresponds to ``gfortran``.
    """
    # NOTE: NumPy may not be installed, but we don't want **this** module to
    #       cause an import failure.
    from numpy.distutils.fcompiler import gnu

    # Only Mac OS X.
    if sys.platform != MAC_OS_X:
        return False

    # Only ``gfortran``.
    if not isinstance(f90_compiler, gnu.Gnu95FCompiler):
        return False

    return True


def is_dual_architecture():
    """Checks if the current Python binary is dual architecture.

    Expected to be used by :func:`patch__compile`, so assumes the
    caller has already checked this is running on OS X.

    This uses ``lipo -info`` to check that the executable is a "fat file"
    with both ``i386`` and ``x86_64`` architectures.

    We use ``lipo -info`` rather than ``file`` because ``lipo`` is
    purpose-built for checking the architecture(s) in a file.

    This property could also be checked by looking for the presence of
    multiple architectures in
    ``distutils.sysconfig.get_config_var('LDFLAGS')``.

    Returns:
        bool: Indicating if the Python binary is dual architecture
        (:data:`True`) or single architecture (:data:`False`).
    """
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

    Expected to be used by :func:`patch__compile`, so assumes the
    caller has already checked this is running on OS X.

    By default, the Homebrew ``gfortran`` **does not** support building dual
    architecture object files. This checks support for this feature by trying
    to build a very simple Fortran 90 program with ``-arch i386 -arch x86_64``.

    Returns:
        bool: Indicates if ``gfortran`` can build dual architecture binaries.
    """
    temporary_directory = tempfile.mkdtemp(suffix='-fortran')
    source_name = os.path.join(temporary_directory, 'bar.f90')
    with open(source_name, 'w') as file_obj:
        file_obj.write(SIMPLE_F90_SOURCE)

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


def patch__compile(f90_compiler):
    """Modify the Fortran compiler to create universal binary object files.

    Assumes the caller is :func:`patch_f90_compiler` (i.e. that
    :func:`is_osx_gfortran` has returned :data:`True`).

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

    f90_compiler._compile = DualArchitectureCompile(f90_compiler)


def patch_f90_compiler(f90_compiler):
    """Patch up ``f90_compiler.library_dirs``.

    On Mac OS X, a Homebrew installed ``gfortran`` has a few quirks that
    must be dealt with.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.
    """
    if not is_osx_gfortran(f90_compiler):
        return

    setup_helpers.patch_library_dirs(f90_compiler)
    patch__compile(f90_compiler)


class DualArchitectureCompile(object):
    """Callable wrapper that over-rides ``_compile``.

    Intended to be used by :func:`patch__compile`.

    Objects of this type are intended to be used to replace / augment to
    ``_compile`` method on a ``Gnu95FCompiler`` (i.e. a Fortran compiler that
    uses ``gfortran``). This is because the Homebrew ``gfortran`` can't build
    fat binaries:

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
        self.compiler_cmds = {
            '.f90': f90_compiler.compiler_f90,
            '.f': f90_compiler.compiler_f77,
        }
        # NOTE: These **should** be cleaned up when this instance goes out
        #       of scope, but that is not particularly simple to do.
        self.i386_dir = self._i386_dir()
        self.x86_64_dir = self._x86_64_dir()

        self.arch_indices = {}  # Populated in ``_verify()``.
        self.arch_values = {}  # Populated in ``_verify()``.
        self._verify()

    def _verify(self):
        """Makes sure the constructor arguments are valid.

        In particular, makes sure that each command in ``compiler_cmds`` has
        exactly one instance of ``-arch``.

        If this succeeds, will update the ``arch_index`` and ``arch_value``
        dictionaries on the instance (for each extension supported in
        ``compiler_cmds``).

        Raises:
            TypeError: If a value in ``compiler_cmds`` is not a ``list``.
            ValueError: If a value in ``compiler_cmds`` doesn't have exactly one
                ``-arch`` segment.
            ValueError: If ``-arch`` is the **last** segment in a value in
                ``compiler_cmds``.
            ValueError: If one of the ``-arch`` values is not ``i386`` or
                ``x86_64``.
        """
        for extension, compiler_cmd in self.compiler_cmds.items():
            if not isinstance(compiler_cmd, list):
                raise TypeError('Expected a list', compiler_cmd, extension)

            if compiler_cmd.count('-arch') != 1:
                raise ValueError(
                    'Did not find exactly one "-arch" in', compiler_cmd,
                    extension)

            arch_index = compiler_cmd.index('-arch') + 1
            if arch_index == len(compiler_cmd):
                raise ValueError(
                    'There is no architecture specified in', compiler_cmd,
                    extension)

            arch_value = compiler_cmd[arch_index]
            if arch_value not in ('i386', 'x86_64'):
                raise ValueError(
                    'Unexpected architecture', arch_value, 'in',
                    compiler_cmd, extension)

            self.arch_indices[extension] = arch_index
            self.arch_values[extension] = arch_value

    def _set_architecture(self, architecture, extension):
        """Set the architecture on the Fortran compiler.

        ``compiler_cmd`` is actually a list (mutable), so we can update it here
        and it will change the architecture that ``f90_compiler`` targets.

        Args:
            architecture (str): One of ``i386`` or ``x86_64``.
            extension (str): The file extension being compiled.
        """
        compiler_cmd = self.compiler_cmds[extension]
        arch_index = self.arch_indices[extension]
        compiler_cmd[arch_index] = architecture

    def _restore_architecture(self, extension):
        """Restore the architecture on the Fortran compiler.

        Resets the ``-arch`` value in ``compiler_cmd`` to its original value.

        Args:
            extension (str): The file extension being compiled.
        """
        compiler_cmd = self.compiler_cmds[extension]
        arch_index = self.arch_indices[extension]
        arch_value = self.arch_values[extension]
        compiler_cmd[arch_index] = arch_value

    @staticmethod
    def _i386_dir():
        """Temporary directory for building objects that target ``i386``.

        If the ``TEMPDIR_I386`` environment variable is set (and that is a
        valid directory), will just return that. Otherwise, will create a
        temporary directory.

        The environment variable is really just intended to be used by Travis
        CI so that the build journal is "consistent".

        Returns:
            str: The directory created (or specified from the environment).
        """
        i386_dir = os.environ.get('TEMPDIR_I386')
        if i386_dir is None or not os.path.isdir(i386_dir):
            i386_dir = tempfile.mkdtemp(suffix='-i386')
        return i386_dir

    @staticmethod
    def _x86_64_dir():
        """Temporary directory for building objects that target ``x86_64``.

        If the ``TEMPDIR_X86_64`` environment variable is set (and that is a
        valid directory), will just return that. Otherwise, will create a
        temporary directory.

        The environment variable is really just intended to be used by Travis
        CI so that the build journal is "consistent".

        Returns:
            str: The directory created (or specified from the environment).
        """
        x86_64_dir = os.environ.get('TEMPDIR_X86_64')
        if x86_64_dir is None or not os.path.isdir(x86_64_dir):
            x86_64_dir = tempfile.mkdtemp(suffix='-x86_64')
        return x86_64_dir

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
        i386_dir = self.i386_dir
        i386_obj = os.path.join(i386_dir, obj_name)
        self._set_architecture('i386', ext)
        self.original_compile(
            i386_obj, src, ext, cc_args, extra_postargs, pp_opts)

        # Create a directory and compile an object targeting x86_64.
        x86_64_dir = self.x86_64_dir
        x86_64_obj = os.path.join(x86_64_dir, obj_name)
        self._set_architecture('x86_64', ext)
        self.original_compile(
            x86_64_obj, src, ext, cc_args, extra_postargs, pp_opts)

        # Restore the compiler back to how it was before we modified it (could
        # be done with a context manager).
        self._restore_architecture(ext)

        # Use ``lipo`` to combine the object files into a universal.
        lipo_cmd = ('lipo', i386_obj, x86_64_obj, '-create', '-output', obj)
        self.f90_compiler.spawn(lipo_cmd)
