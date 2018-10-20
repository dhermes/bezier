#####################
C ABI (``libbezier``)
#####################

The core computational routines used by ``bezier`` are implemented in
`Fortran`_. The subroutines provided there are available in a C ABI
``libbezier``. They can be called from Fortran, C, C++, Cython and any other
language that can invoke a foreign C function (e.g. `Go`_ via ``cgo``).

.. _Fortran: https://en.wikipedia.org/wiki/Fortran
.. _Go: https://golang.org

The modules below correspond to the Fortran modules used to structure the
library. The corresponding header files for the ABI match this module
structure.

.. toctree::
   :titlesonly:

   curve
   curve_intersection
   helpers
   status
   surface
   surface_intersection

************
Installation
************

Currently (as of October 19, 2018) there is no direct support for building
``libbezier`` or for installing the header files.

However, the :doc:`../python-binary-extension` document explains how to locate
a built ``libbezier`` static library and the C headers once the ``bezier``
Python package has been installed.
