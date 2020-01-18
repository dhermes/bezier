############################
Cython ``.pxd`` Declarations
############################

In addition to the header files, several ``cimport``-able ``.pxd``
Cython declaration files are provided:

.. testsetup:: show-pxd

   import os

   import bezier
   import tests.utils


   print_tree = tests.utils.print_tree
   bezier_directory = os.path.dirname(bezier.__file__)

.. doctest:: show-pxd
   :windows-skip:

   >>> bezier_directory
   '.../site-packages/bezier'
   >>> print_tree(bezier_directory, suffix=".pxd")
   bezier/
     _curve.pxd
     _curve_intersection.pxd
     _helpers.pxd
     _status.pxd
     _triangle.pxd
     _triangle_intersection.pxd

For example, ``cimport bezier._curve`` will provide all the functions
in ``bezier/curve.h``.

.. toctree::
   :titlesonly:

   curve <curve>
   curve_intersection <curve_intersection>
   helpers <helpers>
   status <status>
   triangle <triangle>
   triangle_intersection <triangle_intersection>
