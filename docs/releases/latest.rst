Latest Release (``0.8.1.dev1``)
===============================

|pypi| |docs|

ABI Changes
-----------

New Features
~~~~~~~~~~~~

-  Added ``surface.h::compute_area`` helper that can be used to compute the
   area of both a surface and a curved polygon
   (`d4d7249 <https://github.com/dhermes/bezier/commit/d4d7249729dffd4994df1af899084ceb89dde8fc>`__).

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removing getters and setters for parameters used during curve-curve
   intersection
   (`2fda3ae <https://github.com/dhermes/bezier/commit/2fda3aed2818849363c425e3fce70b4bafe7e9ef>`__):

   -  ``curve_intersection.h::set_max_candidates``
   -  ``curve_intersection.h::get_max_candidates``

Python Changes
--------------

New Features
~~~~~~~~~~~~

-  Added implementation for ``Surface.area``
   `property <http://bezier.readthedocs.io/en/0.8.1/reference/bezier.surface.html#bezier.surface.Surface.area>`__
   and ``CurvedPolygon.area``
   `property <http://bezier.readthedocs.io/en/0.8.1/reference/bezier.curved_polygon.html#bezier.curved_polygon.CurvedPolygon.area>`__
   (`eb6077e <https://github.com/dhermes/bezier/commit/eb6077eab4f6ca0d72de6194f1789a2d0eada8b0>`__).

Non-Public API
~~~~~~~~~~~~~~

-  Removing getters and setters for parameters used during curve-curve
   intersection
   (`2fda3ae <https://github.com/dhermes/bezier/commit/2fda3aed2818849363c425e3fce70b4bafe7e9ef>`__):

   -  ``bezier._geometric_intersection.set_max_candidates()``
   -  ``bezier._geometric_intersection.get_max_candidates()``
-  Removing cached values for ``Curve.length``
   `property <http://bezier.readthedocs.io/en/0.8.1/reference/bezier.curve.html#bezier.curve.Curve.length>`__,
   ``Surface.area``
   `property <http://bezier.readthedocs.io/en/0.8.1/reference/bezier.surface.html#bezier.surface.Surface.area>`__
   and ``Surface.is_valid``
   `property <http://bezier.readthedocs.io/en/0.8.1/reference/bezier.surface.html#bezier.surface.Surface.is_valid>`__
   (`34d48d6 <https://github.com/dhermes/bezier/commit/34d48d6900963734d7fb82f13bd3f37416cc6efe>`__).

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.8.1.svg
   :target: https://pypi.org/project/bezier/0.8.1/
   :alt: PyPI link to release 0.8.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.8.1
   :target: https://bezier.readthedocs.io/en/0.8.1/
   :alt: Documentation for release 0.8.1
