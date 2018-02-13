Latest Release (``0.7.1.dev1``)
===============================

|pypi| |docs|

ABI Changes
-----------

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removed ``BAD_TANGENT`` status enum
   (`b89b2b1 <https://github.com/dhermes/bezier/commit/b89b2b1de1726cdc9f508bd761f4c20e7d655321>`__).
   The case where that failure was issued has now been handled as an acceptable
   ``TANGENT_BOTH`` classification for surface-surface intersection points.
   (See the ``classify_intersection()``
   `function <http://bezier.readthedocs.io/en/0.7.1/algorithm-helpers.html#bezier._surface_helpers.classify_intersection>`__
   for an example.)
-  Adding ``BAD_INTERIOR`` status enum
   (`6348dc6 <https://github.com/dhermes/bezier/commit/6348dc63b5d11453fa8312997429448bbdad0a3f>`__).
   (This is a **breaking** change rather than additive because it re-uses
   the enum value of ``5`` previously used by ``BAD_TANGENT``.) This
   value is used by ``interior_combine()`` in the case that the
   curved polygon intersection(s) cannot be determined from the edge-edge
   intersections for a given surface-surface pair. See
   `#101 <https://github.com/dhermes/bezier/issues/101>`__.

Python Changes
--------------

-  (**Bug fix**) The ``0.7.0`` release broke ``Surface.plot()`` and
   ``CurvedPolygon.plot()`` (when the nodes were transposed, the plotting
   helpers were not correctly updated). The ``add_patch()`` helper was
   fixed to account for the changes in data layout
   (`80bfaaa <https://github.com/dhermes/bezier/commit/80bfaaa74219f9053585aa8970131018baa516d1>`__).

Non-Public API
~~~~~~~~~~~~~~

-  Adding a ``TANGENT_BOTH`` classification for surface-surface intersection
   points that are interior to both surfaces at the point of tangency
   (`b89b2b1 <https://github.com/dhermes/bezier/commit/b89b2b1de1726cdc9f508bd761f4c20e7d655321>`__).
   This previously failed with a :exc:`NotImplementedError`.

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.7.1.svg
   :target: https://pypi.org/project/bezier/0.7.1/
   :alt: PyPI link to release 0.7.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.7.1
   :target: https://bezier.readthedocs.io/en/0.7.1/
   :alt: Documentation for release 0.7.1
