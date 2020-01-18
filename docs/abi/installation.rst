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
   >   -DCMAKE_BUILD_TYPE=Release \
   >   -DCMAKE_INSTALL_PREFIX:PATH="${INSTALL_PREFIX}" \
   >   -S "${SRC_DIR}" \
   >   -B "${BUILD_DIR}"
   $ cmake \
   >   --build "${BUILD_DIR}" \
   >   --config Debug \
   >   --target install

Note that this will require having a full checkout of the Fortran source
code (in ``SRC_DIR``).

*******************
Installed Artifacts
*******************

Once installed, the ``libbezier`` library will come with a shared library
for linking and C headers.

Linux
=====

.. code-block:: console

   $ tree "${INSTALL_PREFIX}"
   .../libbezier-release/usr/
   ├── include
   │   ├── bezier
   │   │   ├── curve.h
   │   │   ├── curve_intersection.h
   │   │   ├── helpers.h
   │   │   ├── status.h
   │   │   ├── triangle.h
   │   │   └── triangle_intersection.h
   │   └── bezier.h
   ├── lib
   │   ├── libbezier.so -> libbezier.so.2020
   │   ├── libbezier.so.2020 -> libbezier.so.2020.1.14
   │   └── libbezier.so.2020.1.14
   └── share
       └── bezier
           └── cmake
               ├── BezierConfig-release.cmake
               └── BezierConfig.cmake

   6 directories, 12 files

macOS
=====

.. code-block:: console

   $ tree "${INSTALL_PREFIX}"
   .../libbezier-release/usr/
   ├── include
   │   ├── bezier
   │   │   ├── curve.h
   │   │   ├── curve_intersection.h
   │   │   ├── helpers.h
   │   │   ├── status.h
   │   │   ├── triangle.h
   │   │   └── triangle_intersection.h
   │   └── bezier.h
   ├── lib
   │   ├── libbezier.2020.1.14.dylib
   │   ├── libbezier.2020.dylib -> libbezier.2020.1.14.dylib
   │   └── libbezier.dylib -> libbezier.2020.dylib
   └── share
       └── bezier
           └── cmake
               ├── BezierConfig-release.cmake
               └── BezierConfig.cmake

   6 directories, 12 files

Windows
=======

.. code-block:: console

   $ tree "${INSTALL_PREFIX}"
   .../libbezier-release/usr/
   ├── bin
   │   └── bezier.dll
   ├── include
   │   ├── bezier
   │   │   ├── curve.h
   │   │   ├── curve_intersection.h
   │   │   ├── helpers.h
   │   │   ├── status.h
   │   │   ├── triangle.h
   │   │   └── triangle_intersection.h
   │   └── bezier.h
   ├── lib
   │   └── bezier.lib
   └── share
       └── bezier
           └── cmake
               ├── BezierConfig-release.cmake
               └── BezierConfig.cmake

   7 directories, 11 files

.. _CMake: https://cmake.org/
