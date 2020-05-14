Latest Release (``2020.2.4.dev1``)
==================================

|pypi| |docs|

Python Changes
--------------

Breaking Changes
~~~~~~~~~~~~~~~~

-  Moved non-public ``bezier._algebraic_intersection`` module to
   ``bezier.hazmat.algebraic_intersection``
   (`#216 <https://github.com/dhermes/bezier/pull/216>`__).
-  Moved non-public ``bezier._py_curve_helpers`` module to
   ``bezier.hazmat.curve_helpers``
   (`#218 <https://github.com/dhermes/bezier/pull/218>`__).
-  Moved non-public ``bezier._py_triangle_intersection`` module to
   ``bezier.hazmat.triangle_intersection``
   (`#219 <https://github.com/dhermes/bezier/pull/219>`__).
-  Moved non-public ``bezier._py_triangle_helpers`` module to
   ``bezier.hazmat.triangle_helpers``
   (`#220 <https://github.com/dhermes/bezier/pull/220>`__).
-  Moved non-public ``bezier._py_intersection_helpers`` module to
   ``bezier.hazmat.intersection_helpers``
   (`#222 <https://github.com/dhermes/bezier/pull/222>`__).

Documentation
--------------

-  Removed ``algorithms/algebraic-helpers`` document since the
   ``bezier.hazmat.algebraic_intersection`` module is now fully documented
   (`#216 <https://github.com/dhermes/bezier/pull/216>`__).
-  Updated from ``https://docs.scipy.org/doc/numpy`` to ``https://numpy.org``
   for references to the NumPy documentation
   (`#221 <https://github.com/dhermes/bezier/pull/221>`__).

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/2020.2.4.dev1.svg
   :target: https://pypi.org/project/bezier/2020.2.4.dev1/
   :alt: PyPI link to release 2020.2.4.dev1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=2020.2.4.dev1
   :target: https://bezier.readthedocs.io/en/2020.2.4.dev1/
   :alt: Documentation for release 2020.2.4.dev1
