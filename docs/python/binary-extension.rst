################
Binary Extension
################

.. note::

   This content was last updated August 1, 2023. Much of the content is tested
   automatically to keep from getting stale. If anything does not seem correct
   --- or even if the explanation is insufficient --- please `file an issue`_.

   .. _file an issue: https://github.com/dhermes/bezier/issues/new

The ``bezier`` Python package has optional speedups that wrap the
``libbezier`` library. These are incorporated into the Python interface via
`Cython`_ as a binary extension. See :doc:`../abi/index` for more information
on building and installing ``libbezier``.

.. _Cython: https://cython.readthedocs.io/

***************************
Extra (Binary) Dependencies
***************************

When the ``bezier`` Python package is installed via `pip`_, it will likely be
installed from a `Python wheel`_. The wheels uploaded to PyPI are pre-built,
with the Fortran code compiled by `GNU Fortran`_ (``gfortran``). As a
result, ``libbezier`` will depend on ``libgfortran``. This can be problematic
due to version conflicts, ABI incompatibility, a desire to use a different
Fortran compiler (e.g. Intel's ``ifort``) and a host of other reasons.

There is standard tooling for distributing wheels that address this:

* Linux: `auditwheel`_
* macOS: `delocate`_
* Windows: `delvewheel`_

.. _auditwheel: https://github.com/pypa/auditwheel
.. _delocate: https://github.com/matthew-brett/delocate
.. _delvewheel: https://github.com/adang1345/delvewheel

The tools address it by placing a copy of ``libgfortran`` (and potentially its
dependencies) in the built wheel. This means that libraries that depend on
``libbezier`` may also need to link against these local copies of dependencies.

.. _pip: https://pip.pypa.io
.. _Python wheel: https://wheel.readthedocs.io
.. _GNU Fortran: https://gcc.gnu.org/fortran/

Linux
=====

The command line tool `auditwheel`_ adds a ``bezier.libs`` directory to
``site-packages`` (i.e. it is **next to** ``bezier``) with a modified
``libbezier`` and all of its dependencies (e.g. ``libgfortran``)

.. testsetup:: linux-libs, linux-readelf-py, linux-readelf-lib, macos-dylibs,
               macos-extension, macos-delocated-libgfortran

   import os
   import subprocess

   import bezier
   import tests.utils


   print_tree = tests.utils.print_tree
   base_dir = os.path.abspath(os.path.dirname(bezier.__file__))
   # macOS specific.
   dylibs_directory = os.path.join(base_dir, ".dylibs")
   # Linux specific.
   libs_directory = os.path.abspath(os.path.join(base_dir, os.pardir, "bezier.libs"))


   def invoke_shell(*args, cwd=base_dir, replacements=()):
       print("$ " + " ".join(args))
       # NOTE: We print to the stdout of the doctest, rather than using
       #       ``subprocess.call()`` directly.
       output_bytes = subprocess.check_output(args, cwd=cwd)
       output_str = output_bytes.decode("utf-8")

       for before, after in replacements:
           output_str = output_str.replace(before, after)

       print(output_str, end="")

.. doctest:: linux-libs
   :linux-only:

   >>> libs_directory
   '.../site-packages/bezier.libs'
   >>> print_tree(libs_directory)
   bezier.libs/
     libbezier-631d8eda.so.2023.7.28
     libgfortran-040039e1.so.5.0.0
     libquadmath-96973f99.so.0.0.0

The ``bezier._speedup`` module depends on this local copy of ``libbezier``:

.. testcode:: linux-readelf-py
   :hide:

   invoke_shell("readelf", "-d", "_speedup.cpython-312-x86_64-linux-gnu.so")

.. testoutput:: linux-readelf-py
   :linux-only:
   :pyversion: >= 3.12

   $ readelf -d _speedup.cpython-312-x86_64-linux-gnu.so

   Dynamic section at offset 0x49f000 contains 27 entries:
     Tag        Type                         Name/Value
    0x000000000000000f (RPATH)              Library rpath: [$ORIGIN/../bezier.libs]
    0x0000000000000001 (NEEDED)             Shared library: [libbezier-631d8eda.so.2023.7.28]
    0x0000000000000001 (NEEDED)             Shared library: [libpthread.so.0]
    0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
    0x000000000000000c (INIT)               0x7000
    0x000000000000000d (FINI)               0x809f0
   ...

and the local copy of ``libbezier`` depends on the other dependencies in
``bezier.libs/`` (both directly and indirectly):

.. testcode:: linux-readelf-lib
   :hide:

   invoke_shell("readelf", "-d", "../bezier.libs/libbezier-631d8eda.so.2023.7.28")
   invoke_shell("readelf", "-d", "../bezier.libs/libgfortran-040039e1.so.5.0.0")

.. testoutput:: linux-readelf-lib
   :linux-only:

   $ readelf -d ../bezier.libs/libbezier-631d8eda.so.2023.7.28

   Dynamic section at offset 0x4adc8 contains 29 entries:
     Tag        Type                         Name/Value
    0x0000000000000001 (NEEDED)             Shared library: [libgfortran-040039e1.so.5.0.0]
    0x0000000000000001 (NEEDED)             Shared library: [libm.so.6]
    0x0000000000000001 (NEEDED)             Shared library: [libgcc_s.so.1]
    0x0000000000000001 (NEEDED)             Shared library: [libquadmath-96973f99.so.0.0.0]
    0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
    0x000000000000000e (SONAME)             Library soname: [libbezier-631d8eda.so.2023.7.28]
    0x000000000000000c (INIT)               0x3000
   ...
   $ readelf -d ../bezier.libs/libgfortran-040039e1.so.5.0.0

   Dynamic section at offset 0x275d78 contains 31 entries:
     Tag        Type                         Name/Value
    0x0000000000000001 (NEEDED)             Shared library: [libquadmath-96973f99.so.0.0.0]
    0x0000000000000001 (NEEDED)             Shared library: [libz.so.1]
    0x0000000000000001 (NEEDED)             Shared library: [libm.so.6]
    0x0000000000000001 (NEEDED)             Shared library: [libgcc_s.so.1]
    0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
    0x000000000000000e (SONAME)             Library soname: [libgfortran-040039e1.so.5.0.0]
    0x000000000000000c (INIT)               0x19a88
   ...

.. note::

   The runtime path (``RPATH``) uses ``$ORIGIN`` to specify a path
   relative to the directory where the extension module (``.so`` file) is.

macOS
=====

The command line tool `delocate`_ adds a ``bezier/.dylibs`` directory
with copies of ``libbezier``, ``libgfortran``, ``libquadmath`` and
``libgcc_s``:

.. doctest:: macos-dylibs
   :macos-only:

   >>> dylibs_directory
   '.../site-packages/bezier/.dylibs'
   >>> print_tree(dylibs_directory)
   .dylibs/
     libbezier.2023.7.28.dylib
     libgcc_s.1.1.dylib
     libgfortran.5.dylib
     libquadmath.0.dylib

The ``bezier._speedup`` module depends on the local copy
of ``libbezier``:

.. testcode:: macos-extension
   :hide:

   invoke_shell(
      "otool",
      "-L",
      "_speedup.cpython-312-darwin.so",
      replacements=(("\t", "        "),),
   )

.. testoutput:: macos-extension
   :macos-only:
   :pyversion: >= 3.12

   $ otool -L _speedup.cpython-312-darwin.so
   _speedup.cpython-312-darwin.so:
           @loader_path/.dylibs/libbezier.2023.7.28.dylib (...)
           /usr/lib/libSystem.B.dylib (...)

Though the Python extension module (``.so`` file) only depends on ``libbezier``
it indirectly depends on ``libgfortran``, ``libquadmath`` and ``libgcc_s``:

.. testcode:: macos-delocated-libgfortran
   :hide:

   invoke_shell(
      "otool",
      "-L",
      ".dylibs/libbezier.2023.7.28.dylib",
      replacements=(("\t", "        "),),
   )

.. testoutput:: macos-delocated-libgfortran
   :macos-only:

   $ otool -L .dylibs/libbezier.2023.7.28.dylib
   .dylibs/libbezier.2023.7.28.dylib:
           /DLC/bezier/.dylibs/libbezier.2023.7.28.dylib (...)
           @loader_path/libgfortran.5.dylib (...)
           @loader_path/libquadmath.0.dylib (...)
           /usr/lib/libSystem.B.dylib (...)

.. note::

   To allow the package to be relocatable, the ``libbezier`` dependency is
   relative to the ``@loader_path`` (i.e. the path where the Python extension
   module is loaded) instead of being an absolute path within the file
   system.

   Notice also that ``delocate`` uses the nonexistent root ``/DLC`` for
   the ``install_name`` of ``libbezier`` to avoid accidentally pointing
   to an existing file on the target system.

Windows
=======

The command line tool `delvewheel`_ adds a ``bezier.libs`` directory to
``site-packages`` (i.e. it is **next to** ``bezier``) with a modified
``libbezier`` DLL

.. doctest:: windows-libs
   :windows-only:

   >>> libs_directory
   '...\\site-packages\\bezier.libs'
   >>> print_tree(libs_directory)
   bezier.libs\
     bezier-40ff1ce7372f05ba11436ffbadd11324.dll
     libgcc_s_seh-1-5c71c85c0ca01174917203266ba98140.dll
     libgfortran-5-08073c6868a1df2cbc5609e49cbe3ad8.dll
     libquadmath-0-55d07eaa5b490be06911c864dcae60fd.dll
     libwinpthread-1-737bdf20e708783437e6fdbd7b05edf7.dll

The ``bezier._speedup`` module (``.pyd`` file) depends on this local copy of
``libbezier``:

.. testsetup:: windows-libs, windows-extension, windows-dll

   import distutils.ccompiler
   import os
   import pathlib
   import re
   import subprocess

   import bezier
   import tests.utils


   base_dir = os.path.abspath(os.path.dirname(bezier.__file__))
   site_packages = os.path.abspath(os.path.join(base_dir, os.pardir))
   libs_directory = os.path.join(site_packages, "bezier.libs")
   # Use regex replacement to handle the fact that the ``bezier.dll``
   # file contents are non-deterministic (across time / builds). The
   # MinGW packages **are** deterministic (for a given version) but those
   # may differ across different build machines so we replace them too.
   dll_replacements = (
      ("bezier-[0-9a-f]{32}.dll", "bezier-40ff1ce7372f05ba11436ffbadd11324.dll"),
      ("libgcc_s_seh-1-[0-9a-f]{32}.dll", "libgcc_s_seh-1-5c71c85c0ca01174917203266ba98140.dll"),
      ("libgfortran-5-[0-9a-f]{32}.dll", "libgfortran-5-08073c6868a1df2cbc5609e49cbe3ad8.dll"),
      ("libquadmath-0-[0-9a-f]{32}.dll", "libquadmath-0-55d07eaa5b490be06911c864dcae60fd.dll"),
      ("libwinpthread-1-[0-9a-f]{32}.dll", "libwinpthread-1-737bdf20e708783437e6fdbd7b05edf7.dll"),
   )


   def print_tree(directory):
       return tests.utils.print_tree(directory, replacements=dll_replacements)


   if os.name == "nt":
       c_compiler = distutils.ccompiler.new_compiler()
       assert c_compiler.compiler_type == "msvc"
       c_compiler.initialize()

       dumpbin_exe = os.path.join(
           os.path.dirname(c_compiler.lib), "dumpbin.exe")
       assert os.path.isfile(dumpbin_exe)
   else:
       # This won't matter if not on Windows.
       dumpbin_exe = None


   def replace_dumpbin(value):
       if value == "dumpbin":
           return dumpbin_exe
       else:
           return value


   def _transform_deps_sort_func(value):
       if value.strip().startswith("lib"):
           return 1, value.lower()

       return 2, value.lower()


   def transform_deps(output_str):
       separator_before = "\n  Image has the following dependencies:\n\n"
       separator_after = "\n\n  Summary\n"
       before, partial = output_str.split(separator_before)
       dependency_str, after = partial.split(separator_after)
       dependencies = dependency_str.split("\n")
       dependencies.sort(key=_transform_deps_sort_func)

       modified = (
           before
           + separator_before
           + "\n".join(dependencies)
           + separator_after
           + after
       )
       return modified


   def invoke_shell(*args, cwd=base_dir, transform=None):
       # Replace ``"dumpbin"`` with ``dumpbin_exe``.
       cmd = tuple(map(replace_dumpbin, args))
       # NOTE: We print to the stdout of the doctest, rather than using
       #       ``subprocess.call()`` directly.
       output_bytes = subprocess.check_output(cmd, cwd=cwd)

       output_str = os.linesep.join(
           [
               "> " + " ".join(args),
               output_bytes.decode("utf-8"),
           ]
       )

       for pattern, replacement in dll_replacements:
           output_str = re.sub(pattern, replacement, output_str)

       # Normalize line endings (content is authored with UNIX-style)
       output_str = output_str.replace(os.linesep, "\n")

       if transform is not None:
           output_str = transform(output_str)

       print(output_str, end="")

.. testcode:: windows-extension
   :hide:

   invoke_shell("dumpbin", "/dependents", "_speedup.cp312-win_amd64.pyd")

.. testoutput:: windows-extension
   :windows-only:
   :pyversion: >= 3.12

   > dumpbin /dependents _speedup.cp312-win_amd64.pyd
   Microsoft (R) COFF/PE Dumper Version ...
   Copyright (C) Microsoft Corporation.  All rights reserved.


   Dump of file _speedup.cp312-win_amd64.pyd

   File Type: DLL

     Image has the following dependencies:

       bezier-40ff1ce7372f05ba11436ffbadd11324.dll
       python312.dll
       KERNEL32.dll
       VCRUNTIME140.dll
       api-ms-win-crt-stdio-l1-1-0.dll
       api-ms-win-crt-heap-l1-1-0.dll
       api-ms-win-crt-runtime-l1-1-0.dll
       api-ms-win-crt-math-l1-1-0.dll

     Summary
   ...

and the local copy of ``libbezier`` depends on the other dependencies in
``bezier.libs/`` (both directly and indirectly):

.. testcode:: windows-dll
   :hide:

   site_packages_path = pathlib.Path(site_packages)
   dll_path, = site_packages_path.glob("bezier.libs/bezier-*.dll")
   dll_path = dll_path.relative_to(site_packages_path)
   dll_path = os.path.join(os.pardir, str(dll_path))
   invoke_shell("dumpbin", "/dependents", dll_path, transform=transform_deps)

.. testoutput:: windows-dll
   :windows-only:

   > dumpbin /dependents ..\bezier.libs\bezier-40ff1ce7372f05ba11436ffbadd11324.dll
   Microsoft (R) COFF/PE Dumper Version ...
   Copyright (C) Microsoft Corporation.  All rights reserved.


   Dump of file ..\bezier.libs\bezier-40ff1ce7372f05ba11436ffbadd11324.dll

   File Type: DLL

     Image has the following dependencies:

       libgcc_s_seh-1-5c71c85c0ca01174917203266ba98140.dll
       libgfortran-5-08073c6868a1df2cbc5609e49cbe3ad8.dll
       api-ms-win-crt-environment-l1-1-0.dll
       api-ms-win-crt-heap-l1-1-0.dll
       api-ms-win-crt-math-l1-1-0.dll
       api-ms-win-crt-private-l1-1-0.dll
       api-ms-win-crt-runtime-l1-1-0.dll
       api-ms-win-crt-stdio-l1-1-0.dll
       api-ms-win-crt-string-l1-1-0.dll
       api-ms-win-crt-time-l1-1-0.dll
       KERNEL32.dll

     Summary
   ...

To enable building the Python binary extension, the ``libbezier`` DLL also has
a corresponding `import library`_ --- ``usr/lib/bezier.lib`` --- which is
provided to specify the symbols in the DLL.

.. _import library: https://docs.python.org/3/extending/windows.html#differences-between-unix-and-windows

On Windows, building Python extensions is a bit more constrained. Each
official Python is built with a particular `version of MSVC`_ and
Python extension modules must be built with the same compiler. This
is primarily because the C runtime (provided by Microsoft) **changes** from
Python version to Python version. To see why the same C runtime must be used,
consider the following example. If an extension uses ``malloc`` from
``MSVCRT.dll`` to allocate memory for an object and the Python interpreter
tries to free that memory with ``free`` from ``MSVCR90.dll``, `bad things`_
can happen:

.. _bad things: https://stackoverflow.com/questions/30790494/what-are-the-differences-among-the-ways-to-access-msvcrt-in-python-on-windows#comment49633975_30790494

    Python's linked CRT, which is ``msvcr90.dll`` for Python 2.7,
    ``msvcr100.dll`` for Python 3.4, and several ``api-ms-win-crt`` DLLs
    (forwarded to ``ucrtbase.dll``) for Python 3.5 ... Additionally each CRT
    uses its own heap for malloc and free (wrapping Windows ``HeapAlloc`` and
    ``HeapFree``), so allocating memory with one and freeing with another is
    an error.

This problem has been `largely fixed`_ in newer versions of Python but is
still worth knowing.

Unfortunately, there is no Fortran compiler provided by MSVC. The
`MinGW-w64`_ suite of tools is a port of the GNU Compiler Collection (``gcc``)
for Windows. In particular, MinGW includes ``gfortran``. However, mixing the
two compiler families (MSVC and MinGW) can be problematic because MinGW uses
a fixed version of the C runtime (``MSVCRT.dll``) and this dependency cannot
be easily dropped or changed.

.. note::

   Although ``msvcrt.dll`` is a dependency of ``bezier.dll``, it is not
   a problem. Any values returned from Fortran (as ``intent(out)``) will
   have already been allocated by the caller (e.g. the Python process).
   This won't necessarily be true for generic Fortran subroutines, but
   subroutines marked with ``bind(c)`` (i.e. marked as part of the C ABI
   of ``libbezier``) will not be allowed to use ``allocatable`` or
   `deferred-shape`_ output variables. Any memory allocated in Fortran will be
   isolated within the Fortran code.

   .. _deferred-shape: https://thinkingeek.com/2017/01/14/gfortran-array-descriptor/

.. _version of MSVC: https://matthew-brett.github.io/pydagogue/python_msvc.html
.. _largely fixed: https://stevedower.id.au/blog/building-for-python-3-5-part-two
.. _MinGW-w64: https://mingw-w64.org

Source
======

For code that depends on ``libgfortran``, it may be problematic to **also**
depend on the local copy distributed with the ``bezier`` wheels.

The ``bezier`` Python package can be built from source if it is not feasible to
link with these libraries, if a different Fortran compiler is required or
"just because".

The Python extension module can be built from source via:

.. code-block:: console

   $ # One of
   $ BEZIER_INSTALL_PREFIX=.../usr/ python -m pip wheel .
   $ BEZIER_INSTALL_PREFIX=.../usr/ python -m pip install .
   $ BEZIER_INSTALL_PREFIX=.../usr/ python setup.py build_ext
   $ BEZIER_INSTALL_PREFIX=.../usr/ python setup.py build_ext --inplace
