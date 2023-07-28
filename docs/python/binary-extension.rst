################
Binary Extension
################

.. note::

   This content was last updated May 19, 2020. Much of the content is tested
   automatically to keep from getting stale, but some of the console code
   blocks are not. As a result, this material may be out of date. If anything
   does not seem correct --- or even if the explanation is insufficient ---
   please `file an issue`_.

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

Some of the standard tooling for distributing wheels tries to address this. On
Linux and macOS, the tools address it by placing a copy of ``libgfortran`` (and
potentially its dependencies) in the built wheel. (On Windows, there is no
standard tooling beyond that provided by ``distutils`` and ``setuptools``.)
This means that libraries that depend on ``libbezier`` may also need to link
against these local copies of dependencies.

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
   libs_directory = os.path.abspath(os.path.join(base_dir, "..", "bezier.libs"))


   def invoke_shell(*args):
       print("$ " + " ".join(args))
       # NOTE: We print to the stdout of the doctest, rather than using
       #       ``subprocess.call()`` directly.
       output_bytes = subprocess.check_output(args, cwd=base_dir)
       print(output_bytes.decode("utf-8"), end="")

.. doctest:: linux-libs
   :linux-only:

   >>> libs_directory
   '.../site-packages/bezier.libs'
   >>> print_tree(libs_directory)
   bezier.libs/
     libbezier-631d8eda.so.2023.7.27
     libgfortran-040039e1.so.5.0.0
     libquadmath-96973f99.so.0.0.0

The ``bezier._speedup`` module depends on this local copy of ``libbezier``:

.. testcode:: linux-readelf-py
   :hide:

   invoke_shell("readelf", "-d", "_speedup.cpython-311-x86_64-linux-gnu.so")

.. testoutput:: linux-readelf-py
   :linux-only:

   $ readelf -d _speedup.cpython-311-x86_64-linux-gnu.so

   Dynamic section at offset 0x4a9000 contains 27 entries:
     Tag        Type                         Name/Value
    0x000000000000000f (RPATH)              Library rpath: [$ORIGIN/../bezier.libs]
    0x0000000000000001 (NEEDED)             Shared library: [libbezier-631d8eda.so.2023.7.27]
    0x0000000000000001 (NEEDED)             Shared library: [libpthread.so.0]
    0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
    0x000000000000000c (INIT)               0x6000
    0x000000000000000d (FINI)               0x80c80
   ...

and the local copy of ``libbezier`` depends on the other dependencies in
``bezier.libs/`` (both directly and indirectly):

.. testcode:: linux-readelf-lib
   :hide:

   invoke_shell("readelf", "-d", "../bezier.libs/libbezier-631d8eda.so.2023.7.27")
   invoke_shell("readelf", "-d", "../bezier.libs/libgfortran-040039e1.so.5.0.0")

.. testoutput:: linux-readelf-lib
   :linux-only:

   $ readelf -d ../bezier.libs/libbezier-631d8eda.so.2023.7.27

   Dynamic section at offset 0x4adc8 contains 29 entries:
     Tag        Type                         Name/Value
    0x0000000000000001 (NEEDED)             Shared library: [libgfortran-040039e1.so.5.0.0]
    0x0000000000000001 (NEEDED)             Shared library: [libm.so.6]
    0x0000000000000001 (NEEDED)             Shared library: [libgcc_s.so.1]
    0x0000000000000001 (NEEDED)             Shared library: [libquadmath-96973f99.so.0.0.0]
    0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
    0x000000000000000e (SONAME)             Library soname: [libbezier-631d8eda.so.2023.7.27]
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

.. _auditwheel: https://github.com/pypa/auditwheel

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
     libbezier.2023.7.27.dylib
     libgcc_s.1.1.dylib
     libgfortran.5.dylib
     libquadmath.0.dylib

The ``bezier._speedup`` module depends on the local copy
of ``libbezier``:

.. testcode:: macos-extension
   :hide:

   invoke_shell("otool", "-L", "_speedup.cpython-311-darwin.so")

.. testoutput:: macos-extension
   :options: +NORMALIZE_WHITESPACE
   :macos-only:
   :pyversion: >= 3.11

   $ otool -L _speedup.cpython-311-darwin.so
   _speedup.cpython-311-darwin.so:
           @loader_path/.dylibs/libbezier.2023.7.27.dylib (...)
           /usr/lib/libSystem.B.dylib (...)

Though the Python extension module (``.so`` file) only depends on ``libbezier``
it indirectly depends on ``libgfortran``, ``libquadmath`` and ``libgcc_s``:

.. testcode:: macos-delocated-libgfortran
   :hide:

   invoke_shell("otool", "-L", ".dylibs/libbezier.2023.7.27.dylib")

.. testoutput:: macos-delocated-libgfortran
   :options: +NORMALIZE_WHITESPACE
   :macos-only:

   $ otool -L .dylibs/libbezier.2023.7.27.dylib
   .dylibs/libbezier.2023.7.27.dylib:
       /DLC/bezier/.dylibs/libbezier.2023.7.27.dylib (...)
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

.. _delocate: https://github.com/matthew-brett/delocate

Windows
=======

A single Windows shared library (DLL) is provided: ``bezier.dll``.
The Python extension module (``.pyd`` file) depends directly on this library:

.. testsetup:: windows-extension, windows-dll

   import distutils.ccompiler
   import os
   import subprocess

   import bezier


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

   bezier_directory = os.path.dirname(bezier.__file__)


   def replace_dumpbin(value):
       if value == "dumpbin":
           return dumpbin_exe
       else:
           return value


   def invoke_shell(*args):
       print("> " + " ".join(args))
       # Replace ``"dumpbin"`` with ``dumpbin_exe``.
       cmd = tuple(map(replace_dumpbin, args))
       # NOTE: We print to the stdout of the doctest, rather than using
       #       ``subprocess.call()`` directly.
       output_bytes = subprocess.check_output(cmd, cwd=bezier_directory)
       print(output_bytes.decode("utf-8"), end="")

.. testcode:: windows-extension
   :hide:

   invoke_shell("dumpbin", "/dependents", "_speedup.cp311-win_amd64.pyd")

.. testoutput:: windows-extension
   :options: +NORMALIZE_WHITESPACE
   :windows-only:
   :pyversion: >= 3.11

   > dumpbin /dependents _speedup.cp311-win_amd64.pyd
   Microsoft (R) COFF/PE Dumper Version ...
   Copyright (C) Microsoft Corporation.  All rights reserved.


   Dump of file _speedup.cp311-win_amd64.pyd

   File Type: DLL

     Image has the following dependencies:

       bezier-e5dbb97a.dll
       python311.dll
       KERNEL32.dll
       VCRUNTIME140.dll
       api-ms-win-crt-stdio-l1-1-0.dll
       api-ms-win-crt-heap-l1-1-0.dll
       api-ms-win-crt-runtime-l1-1-0.dll
   ...

For built wheels, the dependency will be renamed from ``bezier.dll`` to a
unique name containing the first 8 characters of the SHA256 hash of the DLL
file (to avoid a name collision) and placed in a directory within the
``bezier`` package: for example ``extra-dll/bezier-e5dbb97a.dll``.
TODO

The ``libbezier`` DLL has **no external dependencies**, but does have
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

A Windows shared library (DLL) can be created after compiling
each of the Fortran submodules:

.. code-block:: console

   $ gfortran \
   >   -shared \
   >   -o bezier.dll \
   >   ${OBJ_FILES} \
   >   -Wl,--output-def,bezier.def

.. note::

   Invoking ``gfortran`` **can** be done from the Windows command prompt, but
   it is easier to do from a shell that explicitly supports MinGW, such as
   MSYS2.

By default, the created shared library will depend on ``gcc`` libraries
provided by MinGW:

.. code-block:: rest

   > dumpbin /dependents ...\bezier.dll
   ...
     Image has the following dependencies:

       KERNEL32.dll
       msvcrt.dll
       libgcc_s_seh-1.dll
       libgfortran-3.dll

Unlike Linux and macOS, on Windows relocating and copying any dependencies
on MinGW (at either compile, link or run time) is explicitly avoided. By adding
the ``-static`` flag

.. code-block:: console
   :emphasize-lines: 2

   $ gfortran \
   >   -static \
   >   -shared \
   >   -o bezier.dll \
   >   ${OBJ_FILES} \
   >   -Wl,--output-def,bezier.def

all the symbols used from ``libgfortran`` or ``libgcc_s`` are statically
included and the resulting shared library ``bezier.dll`` has no dependency
on MinGW:

.. testcode:: windows-dll
   :hide:

   invoke_shell("dumpbin", "/dependents", "extra-dll\\bezier-e5dbb97a.dll")

.. testoutput:: windows-dll
   :options: +NORMALIZE_WHITESPACE
   :windows-only:

   > dumpbin /dependents extra-dll\bezier-e5dbb97a.dll
   Microsoft (R) COFF/PE Dumper Version ...
   Copyright (C) Microsoft Corporation.  All rights reserved.


   Dump of file extra-dll\bezier-e5dbb97a.dll

   File Type: DLL

     Image has the following dependencies:

       KERNEL32.dll
       msvcrt.dll

   Summary
   ...

.. note::

   Although ``msvcrt.dll`` is a dependency of ``bezier.dll``, it is not
   a problem. Any values returned from Fortran (as ``intent(out)``) will
   have already been allocated by the caller (e.g. the Python process).
   This won't necessarily be true for generic Fortran subroutines, but
   subroutines marked with ``bind(c)`` (i.e. marked as part of the C ABI
   of ``libbezier``) will not be allowed to use ``allocatable`` or
   `deferred-shape`_ output variables. Any memory allocated in Fortran will be
   isolated within the Fortran code.

   .. _deferred-shape: http://thinkingeek.com/2017/01/14/gfortran-array-descriptor/

   However, the dependency on ``msvcrt.dll`` can still be avoided if desired.
   The MinGW ``gfortran`` default "specs file" can be captured:

   .. code-block:: console

      $ gfortran -dumpspecs > ${SPECS_FILENAME}

   and modified to replace instances of ``-lmsvcrt`` with a substitute, e.g.
   ``-lmsvcr90``. Then ``gfortran`` can be invoked with the flag
   ``-specs=${SPECS_FILENAME}`` to use the custom spec. (Some
   `other dependencies`_ may also indirectly depend on ``msvcrt.dll``,
   such as ``-lmoldname``. `Removing dependencies`_ is not an easy process.)

   .. _other dependencies: https://www.spiria.com/en/blog/desktop-software/building-mingw-w64-toolchain-links-specific-visual-studio-runtime-library
   .. _Removing dependencies: http://www.pygame.org/wiki/PreparingMinGW

From there, an `import library`_ must be created

.. code-block:: rest

   > lib /def:.\bezier.def /out:.\lib\bezier.lib /machine:${ARCH}

.. note::

   ``lib.exe`` is used from the same version of MSVC that compiled the
   target Python. Luckily ``distutils`` enables this without difficulty.

.. _version of MSVC: http://matthew-brett.github.io/pydagogue/python_msvc.html
.. _largely fixed: http://stevedower.id.au/blog/building-for-python-3-5-part-two/
.. _MinGW-w64: http://mingw-w64.org

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
