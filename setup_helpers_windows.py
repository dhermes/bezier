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

"""Helpers for ``setup.py`` specific to Windows."""


import os
import shutil
import sys

import setup_helpers


LIB_DIR = os.path.join('bezier', 'lib')
DLL_DIR = os.path.join('bezier', 'extra-dll')
DLL_NAME = 'libbezier.dll'
DEF_NAME = 'libbezier.def'
LIB_NAME = 'bezier.lib'
# See: https://docs.python.org/3/library/platform.html#cross-platform
if sys.maxsize == 2**63 - 1:
    MACHINE_TYPE = '/machine:x64'
elif sys.maxsize == 2**31 - 1:
    MACHINE_TYPE = '/machine:x86'
else:  # pragma: NO COVER
    raise ImportError('Unexpected maxsize', sys.maxsize)


def _ensure_exists(*dir_names):
    for dir_name in dir_names:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)


def make_static_lib(build_ext_cmd, obj_files):
    f90_compiler = build_ext_cmd.F90_COMPILER
    c_compiler = f90_compiler.c_compiler
    gfortran_exe =  f90_compiler.compiler_f90[0]

    static_lib_dir = os.path.join(build_ext_cmd.build_lib, LIB_DIR)
    extra_dll_dir = os.path.join(build_ext_cmd.build_lib, DLL_DIR)
    _ensure_exists(static_lib_dir, extra_dll_dir, build_ext_cmd.build_temp)

    # Add build lib dir to extensions (for linking).
    for extension in build_ext_cmd.extensions:
        library_dirs = extension.library_dirs
        if static_lib_dir not in library_dirs:
            library_dirs.append(static_lib_dir)

    # NOTE: Would prefer to use a `tuple` but command must be
    #       mutable on Windows.
    temp_dll_filename = os.path.join(build_ext_cmd.build_temp, DLL_NAME)
    cmd = [
        gfortran_exe,
        '-static',
         '-shared',
         '-o',
         temp_dll_filename,
    ]
    cmd.extend(obj_files)
    def_filename = os.path.join(build_ext_cmd.build_temp, DEF_NAME)
    cmd.append('-Wl,--output-def,' + def_filename)
    f90_compiler.spawn(cmd)

    # NOTE: This assumes, but does not check that ``c_compiler.initialized``
    #       is True.
    temp_lib_filename = os.path.join(build_ext_cmd.build_temp, LIB_NAME)
    cmd = [
        c_compiler.lib,
        '/def:' + def_filename,
        '/out:' + temp_lib_filename,
        MACHINE_TYPE,
    ]
    f90_compiler.spawn(cmd)

    # Move files for distribution from build temp. dir to build lib dir.
    dll_filename = os.path.join(extra_dll_dir, DLL_NAME)
    shutil.move(temp_dll_filename, dll_filename)
    lib_filename = os.path.join(static_lib_dir, LIB_NAME)
    shutil.move(temp_lib_filename, lib_filename)


def run_cleanup(build_ext_cmd):
    """Cleanup after ``BuildFortranThenExt.run``.

    For in-place builds, moves the built shared library into the source
    directory.
    """
    if not build_ext_cmd.inplace:
        return

    bezier_dir = os.path.join('src', 'bezier')
    shutil.move(
        os.path.join(build_ext_cmd.build_lib, LIB_DIR),
        bezier_dir,
    )
    shutil.move(
        os.path.join(build_ext_cmd.build_lib, DLL_DIR),
        bezier_dir,
    )


def patch_f90_compiler(f90_compiler):
    """Patch up ``f90_compiler.library_dirs``.

    Updates flags in ``gfortran`` and ignores other compilers. The only
    modification is the removal of ``-fPIC`` since it is not used on Windows
    and the build flags turn warnings into errors.

    Args:
        f90_compiler (numpy.distutils.fcompiler.FCompiler): A Fortran compiler
            instance.
    """
    # NOTE: NumPy may not be installed, but we don't want **this** module to
    #       cause an import failure.
    from numpy.distutils.fcompiler import gnu

    # Only Windows.
    if os.name != 'nt':
        return

    # Only ``gfortran``.
    if not isinstance(f90_compiler, gnu.Gnu95FCompiler):
        return

    f90_compiler.compiler_f90[:] = [
        value for value in f90_compiler.compiler_f90
        if value != setup_helpers.FPIC
    ]

    c_compiler = f90_compiler.c_compiler
    if c_compiler.compiler_type != 'msvc':
        raise NotImplementedError(
            'MSVC is the only supported C compiler on Windows.')


def patch_cmd(cmd_class):
    # Only Windows.
    if os.name != 'nt':
        return

    cmd_class.CUSTOM_STATIC_LIB = make_static_lib
    cmd_class.USE_SHARED_LIBRARY = True
    cmd_class.CLEANUP = run_cleanup
