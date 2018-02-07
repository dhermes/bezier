Latest Release (``0.7.0``)
==========================

|pypi| |docs|

Robustness
----------

-  Geometric curve-curve intersection has better handling for cases when
   the number of intersection candidates grows large (``MAX_CANDIDES == 64``):

   -  First tries to reduce the number of candidates by checking if the
      **actual** convex hulls of each segment in a candidate pair intersect.
      This is a much "finer" check than using the "blunt" bounding box check.
   -  If the convex hull refinement fails, checks if the curves are coincident,
      i.e. different segments / parameterizations along the same algebraic
      curve. This is done by using the ``Curve.locate()``
      `function <https://bezier.readthedocs.io/en/0.7.0/reference/bezier.curve.html#bezier.curve.Curve.locate>`__
      to try to project each of the four endpoints onto the other curve and
      then re-parameterizing each curve onto a common interval.

Data Structures
---------------

-  Storing ``xy``-points as columns (rather than rows). This was a
   **very** large and breaking change, started in
   `b44af8c <https://github.com/dhermes/bezier/commit/b44af8c3d590add947f905f2bc016af7272fc8e0>`__.
   See `#51 <https://github.com/dhermes/bezier/issues/51>`__ for more
   information.

Python Changes
--------------

Non-Public API
~~~~~~~~~~~~~~

-  Requiring **contiguous** 1D arrays for Cython functions
   (`9ede37d <https://github.com/dhermes/bezier/commit/9ede37dcbb7eda9899a02675939eb4dd66af8e8c>`__).

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.7.0.svg
   :target: https://pypi.org/project/bezier/0.7.0/
   :alt: PyPI link to release 0.7.0
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.7.0
   :target: https://bezier.readthedocs.io/en/0.7.0/
   :alt: Documentation for release 0.7.0
