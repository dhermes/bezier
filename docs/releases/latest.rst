Latest Release (``2021.2.13.dev1``)
===================================

|pypi| |docs|

Robustness
----------

-  Fall back to de Casteljau algorithm for evaluation on high degree curves
   (`#264 <https://github.com/dhermes/bezier/pull/264>`__). Fixed
   `#263 <https://github.com/dhermes/bezier/issues/263>`__. The more efficient
   ``evaluate_multi_vs()``
   `function <https://bezier.readthedocs.io/en/2021.2.13.dev1/python/reference/bezier.hazmat.curve_helpers.html#bezier.hazmat.curve_helpers.evaluate_multi_vs>`__
   relies on the ability to compute binomial coefficients ``(n C k)`` but
   after ``(55 C 26)`` these values cannot be computed in floating
   point without roundoff and for very large degrees these values overflow.

Python Changes
--------------

Bug Fixes
~~~~~~~~~

-  Allow the ``extra-dll`` directory to be absent on Windows. For binary wheel
   installs, this directory contains the ``libbezier`` DLL (e.g.
   ``extra-dll\bezier-2a44d276.dll`` was in the most recent release for 64-bit
   Python 3.9). For pure Python installs, the ``extra-dll`` directory will
   be absent.
   (`#255 <https://github.com/dhermes/bezier/pull/255>`__). Fixed
   `#254 <https://github.com/dhermes/bezier/issues/254>`__.
-  Bug fix for ``clip_range()``
   (`doc <https://bezier.readthedocs.io/en/2021.2.13.dev1/python/reference/bezier.hazmat.clipping.html#bezier.hazmat.clipping.clip_range>`__).
   Intersections with the ``t=0`` / ``t=1`` sides of the fat line were not
   taken into account **and** a value for ``s_max`` could be accidentally left
   unset if the intersections occurred in an order where the ``s``-value
   decreased throughout
   (`#259 <https://github.com/dhermes/bezier/pull/259>`__).

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removing ``Surface`` alias for the ``Triangle``
   `type <https://bezier.readthedocs.io/en/2021.2.13.dev1/python/reference/bezier.triangle.html#bezier.triangle.Triangle>`__
   (`#252 <https://github.com/dhermes/bezier/pull/252>`__). The ``Surface``
   type was deprecated (and converted to an alias) in the ``2020.1.14``
   release.

Additive Changes
~~~~~~~~~~~~~~~~

-  Renaming all "private" ``_verify`` args to ``verify``
   (`#251 <https://github.com/dhermes/bezier/pull/251>`__). For example, in
   ``Curve.intersect()``
   (`doc <https://bezier.readthedocs.io/en/2021.2.13.dev1/python/reference/bezier.curve.html#bezier.curve.Curve.intersect>`__)
-  Adding provisional support for self-intersection checks in planar curves
   (`#265 <https://github.com/dhermes/bezier/pull/265>`__,
   `#267 <https://github.com/dhermes/bezier/pull/267>`__). Fixed
   `#165 <https://github.com/dhermes/bezier/issues/165>`__ and
   `#171 <https://github.com/dhermes/bezier/issues/171>`__.
   Supported via ``Curve.self_intersections()``
   `method <https://bezier.readthedocs.io/en/2021.2.13.dev1/python/reference/bezier.curve.html#bezier.curve.Curve.self_intersections>`__.

Documentation
-------------

-  Making all docs pages import external packages at least once
   (`#257 <https://github.com/dhermes/bezier/pull/257>`__). Fixed
   `#210 <https://github.com/dhermes/bezier/issues/210>`__.

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/2021.2.13.dev1.svg
   :target: https://pypi.org/project/bezier/2021.2.13.dev1/
   :alt: PyPI link to release 2021.2.13.dev1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=2021.2.13.dev1
   :target: https://bezier.readthedocs.io/en/2021.2.13.dev1/
   :alt: Documentation for release 2021.2.13.dev1
