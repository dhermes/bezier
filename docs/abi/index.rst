#####################
C ABI (``libbezier``)
#####################

The core computational routines used by the ``bezier`` Python package are
implemented in `Fortran`_. The subroutines provided there are available in a
C ABI ``libbezier``. They can be called from Fortran, C, C++, Cython and any
other language that can invoke a foreign C function (e.g. `Go`_ via ``cgo``).

.. _Fortran: https://en.wikipedia.org/wiki/Fortran
.. _Go: https://golang.org

The modules below correspond to the Fortran modules used to structure the
library. The corresponding header files for the ABI match this module
structure.

.. toctree::
   :titlesonly:

   installation
   curve
   curve_intersection
   helpers
   status
   triangle
   triangle_intersection
