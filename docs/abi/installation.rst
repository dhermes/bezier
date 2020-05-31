############
Installation
############

To build the Fortran shared library directly, use `CMake`_ version
``3.1`` or later. By `default`_, the ``CMAKE_INSTALL_PREFIX`` will install
the package into ``/usr/local/``. For a more temporary install location,
``-DCMAKE_INSTALL_PREFIX:PATH`` can be provided at build time.

.. _default: https://cmake.org/cmake/help/v3.16/variable/CMAKE_INSTALL_PREFIX.html

.. code-block:: console

   $ SRC_DIR="src/fortran/"
   $ BUILD_DIR=".../libbezier-release/build"
   $ INSTALL_PREFIX=".../libbezier-release/usr"
   $ mkdir -p "${BUILD_DIR}"
   $ cmake \
   >     -DCMAKE_BUILD_TYPE=Release \
   >     -DCMAKE_INSTALL_PREFIX:PATH="${INSTALL_PREFIX}" \
   >     -S "${SRC_DIR}" \
   >     -B "${BUILD_DIR}"
   $ cmake \
   >     --build "${BUILD_DIR}" \
   >     --config Release \
   >     --target install

Note that this will require having a full checkout of the Fortran source
code (in ``SRC_DIR``).

*******************
Installed Artifacts
*******************

Once installed, the ``libbezier`` library will come with a shared library
for linking and C headers.

Linux
=====

.. testsetup:: libbezier-linux, libbezier-macos, libbezier-windows

   import os

   import tests.utils


   install_prefix = os.environ["BEZIER_INSTALL_PREFIX"]
   print_tree = tests.utils.print_tree

.. doctest:: libbezier-linux
   :options: +NORMALIZE_WHITESPACE
   :linux-only:

   >>> print_tree(install_prefix)
   usr/
     include/
       bezier/
         curve.h
         curve.hpp
         curve_intersection.h
         helpers.h
         helpers.hpp
         status.h
         triangle.h
         triangle_intersection.h
       bezier.h
       bezier.hpp
     lib/
       libbezier.so -> libbezier.so.2020
       libbezier.so.2020 -> libbezier.so.2020.5.19
       libbezier.so.2020.5.19
     share/
       bezier/
         cmake/
           BezierConfig-release.cmake
           BezierConfig.cmake

macOS
=====

.. doctest:: libbezier-macos
   :options: +NORMALIZE_WHITESPACE
   :macos-only:

   >>> print_tree(install_prefix)
   usr/
     include/
       bezier/
         curve.h
         curve.hpp
         curve_intersection.h
         helpers.h
         helpers.hpp
         status.h
         triangle.h
         triangle_intersection.h
       bezier.h
       bezier.hpp
     lib/
       libbezier.2020.5.19.dylib
       libbezier.2020.dylib -> libbezier.2020.5.19.dylib
       libbezier.dylib -> libbezier.2020.dylib
     share/
       bezier/
         cmake/
           BezierConfig-release.cmake
           BezierConfig.cmake

Windows
=======

.. doctest:: libbezier-windows
   :options: +NORMALIZE_WHITESPACE
   :windows-only:

   >>> print_tree(install_prefix)
   usr\
     bin\
       bezier.dll
     include\
       bezier\
         curve.h
         curve.hpp
         curve_intersection.h
         helpers.h
         helpers.hpp
         status.h
         triangle.h
         triangle_intersection.h
       bezier.h
       bezier.hpp
     lib\
       bezier.lib
     share\
       bezier\
         cmake\
           BezierConfig-release.cmake
           BezierConfig.cmake

.. _CMake: https://cmake.org/
