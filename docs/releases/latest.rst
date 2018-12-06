Latest Release (``0.9.0``)
==========================

|pypi| |docs|

Documentation
-------------

-  Documenting the C ABI ``libbezier``
   (`4608364 <https://github.com/dhermes/bezier/commit/4608364e9c0a2b3888f7f661e629fceda9d9a431>`__).
   Fixed `#63 <https://github.com/dhermes/bezier/issues/63>`__. This
   `documentation <http://bezier.readthedocs.io/en/0.9.0/abi/index.html>`__
   contains a page for each "module" which corresponds to the underlying
   Fortran module. Each module documents the routines in the corresponding
   header file, e.g. the
   `surface <http://bezier.readthedocs.io/en/0.9.0/abi/surface.html>`__
   document corresponds to the ``bezier/surface.h`` header. Fully working
   C examples have been added for each routine in ``bezier/curve.h`` and for
   the enum in ``bezier/status.h``.
-  Adding section about environment variables to
   `development <http://bezier.readthedocs.io/en/0.9.0/development.html>`__
   document
   (`5186e24 <https://github.com/dhermes/bezier/commit/5186e24a7c7eab5d65ac41ba53e3826b693fc86f>`__).
   Fixed `#78 <https://github.com/dhermes/bezier/issues/78>`__.
-  Remove dependency on ``rawgit.com``
   (`04d0f8d <https://github.com/dhermes/bezier/commit/04d0f8d3155a22c5a048f52f75a3c6ffcc7eba69>`__).
   The website is being turned down. Fixed
   `#130 <https://github.com/dhermes/bezier/issues/130>`__.
-  Renaming the "Native Libraries" document as "Binary Extension"
   (`f99db20 <https://github.com/dhermes/bezier/commit/f99db20312bb4ba7e5943195020a8ced4be9457b>`__).
   In the process, changed most references to the "native" Python extension to
   instead call it a "binary" extension.
-  Added a "Cython ``.pxd``
   `Declarations <http://bezier.readthedocs.io/en/0.9.0/python/pxd/index.html>`__"
   document
   (`f99db20 <https://github.com/dhermes/bezier/commit/f99db20312bb4ba7e5943195020a8ced4be9457b>`__).
   Fixed `#122 <https://github.com/dhermes/bezier/issues/122>`__.
-  Moving all Python specific documentation under a specific URL path
   (`3db483b <https://github.com/dhermes/bezier/commit/3db483b58e2c5dd0f618c15fc01710ec6b1a2907>`__).
   In particular, moving

   -  ``/reference/...`` to ``/python/reference/...``
   -  ``/python-binary-extension.html`` to ``/python/binary-extension.html``
   -  ``/pxd/...`` to ``/python/pxd/...``.

-  Moving all algorithm specific documentation under a specific URL path
   (`6e9c825 <https://github.com/dhermes/bezier/commit/6e9c82501a222c95c616658e6e5e7bc00c9f4288>`__).
   In particular, moving

   -  ``/algorithm-helpers.html`` to ``/algorithms/helpers.html``
   -  ``/curve-curve-intersection.html`` to
      ``/algorithms/curve-curve-intersection.html``

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

-  Removing ``dimension`` from ``curve.h::get_curvature``
   (`1e39c0c <https://github.com/dhermes/bezier/commit/1e39c0ce0502919d83a81902c8d9affdb6c6b892>`__).

Python Changes
--------------

New Features
~~~~~~~~~~~~

-  Added implementation for ``Surface.area``
   `property <http://bezier.readthedocs.io/en/0.9.0/python/reference/bezier.surface.html#bezier.surface.Surface.area>`__
   and ``CurvedPolygon.area``
   `property <http://bezier.readthedocs.io/en/0.9.0/python/reference/bezier.curved_polygon.html#bezier.curved_polygon.CurvedPolygon.area>`__
   (`eb6077e <https://github.com/dhermes/bezier/commit/eb6077eab4f6ca0d72de6194f1789a2d0eada8b0>`__).

Non-Public API
~~~~~~~~~~~~~~

-  Removing getters and setters for parameters used during curve-curve
   intersection
   (`2fda3ae <https://github.com/dhermes/bezier/commit/2fda3aed2818849363c425e3fce70b4bafe7e9ef>`__):

   -  ``bezier._geometric_intersection.set_max_candidates()``
   -  ``bezier._geometric_intersection.get_max_candidates()``
-  Removing cached values for ``Curve.length``
   `property <http://bezier.readthedocs.io/en/0.9.0/python/reference/bezier.curve.html#bezier.curve.Curve.length>`__,
   ``Surface.area``
   `property <http://bezier.readthedocs.io/en/0.9.0/python/reference/bezier.surface.html#bezier.surface.Surface.area>`__
   and ``Surface.is_valid``
   `property <http://bezier.readthedocs.io/en/0.9.0/python/reference/bezier.surface.html#bezier.surface.Surface.is_valid>`__
   (`34d48d6 <https://github.com/dhermes/bezier/commit/34d48d6900963734d7fb82f13bd3f37416cc6efe>`__).

Build
~~~~~

-  Renaming ``libbezier.dll`` shared library to ``bezier.dll`` on Windows
   (`d17a9bc <https://github.com/dhermes/bezier/commit/d17a9bcee194edc9f103734e35023d178ed8923b>`__).
   This follows the correct convention on Windows.
-  Adding Python 3.7 support and making it the default version used for testing
   (`e368e9f <https://github.com/dhermes/bezier/commit/e368e9fd9ab31cfd818fcb9e777dff6dcbd3a7e6>`__).
-  Dropping support for Python 3.5
   (`f99db20 <https://github.com/dhermes/bezier/commit/f99db20312bb4ba7e5943195020a8ced4be9457b>`__).
-  Adding back ``-march=native`` for non-wheel builds
   (`1566019 <https://github.com/dhermes/bezier/commit/1566019635b8ffb8a2e4725a2d51830351e03fa5>`__).
   This way, when installing from source (either via a local checkout or from
   the source distribution on PyPI) the most optimal machine instructions will
   be produced. Fixed `#99 <https://github.com/dhermes/bezier/issues/99>`__.
-  Removing all traces of 32-bit support for macOS
   (`d7620ad <https://github.com/dhermes/bezier/commit/d7620adb862ed6f9be9d2615916f789c3c24c52f>`__).
   This was driven by a
   `decision <https://github.com/numpy/numpy/issues/11625>`__ from the NumPy
   maintainers.

Miscellany
~~~~~~~~~~

-  Adopted ``black``
   `code formatter <https://black.readthedocs.io/en/stable/>`__
   (`f21b52d <https://github.com/dhermes/bezier/commit/f21b52d562daf6c86ddaba326aeee8362361e20f>`__).
-  Adding project URLs and keywords for PyPI
   (`cfb070d <https://github.com/dhermes/bezier/commit/cfb070d651fba4e7df06216a159f623d57036f02>`__).
-  Added 20 new surface-surface functional tests
   (`9fd9c1e <https://github.com/dhermes/bezier/commit/9fd9c1e26138034539e91aed04c97ec497a9e4b2>`__).
   See `#121 <https://github.com/dhermes/bezier/issues/121>`__ for more
   information.
-  Removed time and memory benchmarks due to flakiness and lack of an
   environment that could be used for benchmarking
   (`6a30dc2 <https://github.com/dhermes/bezier/commit/6a30dc22abefe7f7573048659b00fbcd968b8ccc>`__).
   See `#125 <https://github.com/dhermes/bezier/issues/125>`__ to follow
   discussion on re-enabling such benchmarks.
-  Using ``DEBUG=True`` environment variable when running unit tests and
   other related tests
   (`d84dffb <https://github.com/dhermes/bezier/commit/d84dffb9d0e6fe1ee653e01cb9d4297f83aa11e0>`__).

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.9.0.svg
   :target: https://pypi.org/project/bezier/0.9.0/
   :alt: PyPI link to release 0.9.0
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.9.0
   :target: https://bezier.readthedocs.io/en/0.9.0/
   :alt: Documentation for release 0.9.0
