Latest Release (``0.6.1``)
==========================

|pypi| |docs|

Python Changes
--------------

Documentation
~~~~~~~~~~~~~

-  Noting that ``Surface.intersect()`` can return a list of either
   ``CurvedPolygon`` or ``Surface`` instances
   (`16e77d7 <https://github.com/dhermes/bezier/commit/16e77d74c526a216c0c2a74d4536cd1d9f93bcff>`__).

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removing ``IntersectionClassification`` enum from ``_status.pxd``
   (`4da969e <https://github.com/dhermes/bezier/commit/4da969e65cec37ca5c0a56e956e7a1546be24236>`__).

Non-Public API
~~~~~~~~~~~~~~

-  Adding getters and setters for parameters used during curve-curve
   intersection
   (`ef4ebc0 <https://github.com/dhermes/bezier/commit/ef4ebc0654d863610df982f218449b27bd135afc>`__):

   -  ``bezier._geometric_intersection.set_max_candidates()``
   -  ``bezier._geometric_intersection.get_max_candidates()``
   -  ``bezier._geometric_intersection.set_similar_ulps()``
   -  ``bezier._geometric_intersection.get_similar_ulps()``

ABI Changes
-----------

Surface Changes
~~~~~~~~~~~~~~~

-  Switching from ``int`` to an actual enum for relevant functions with
   output values that are enums:

   -  In ``surface_intersection.h::surface_intersections``, ``contained``
      is now a ``SurfaceContained``
      (`0a9c0c3 <https://github.com/dhermes/bezier/commit/0a9c0c3736e95deedeecb8d10284c92ebd39469d>`__)
      and ``status`` is now a ``Status``
      (`c356c32 <https://github.com/dhermes/bezier/commit/c356c32b33781b03785b8868f59efd6ad3076a51>`__)
   -  In ``curve_intersection.h::bbox_intersect``, ``enum_`` is now a
      ``BoxIntersectionType``
      (`ef856af <https://github.com/dhermes/bezier/commit/ef856aff4e87ab0620d1ce28e7fdbd3395c8ec38>`__)
   -  In ``curve_intersection.h::curve_intersections``, ``status`` is now a
      ``Status``
      (`ef856af <https://github.com/dhermes/bezier/commit/ef856aff4e87ab0620d1ce28e7fdbd3395c8ec38>`__)

-  Adding getters and setters for parameters used during curve-curve
   intersection
   (`ef4ebc0 <https://github.com/dhermes/bezier/commit/ef4ebc0654d863610df982f218449b27bd135afc>`__):

   -  ``curve_intersection.h::set_max_candidates``
   -  ``curve_intersection.h::get_max_candidates``
   -  ``curve_intersection.h::set_similar_ulps``
   -  ``curve_intersection.h::get_similar_ulps``

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removing inputs ``curve_start / curve_end`` and outputs
   ``true_start / true_end`` in ``curve.h::specialize_curve``
   (`959c547 <https://github.com/dhermes/bezier/commit/959c5473e97e80b1b4e4fd0109f7e79cf1dc36eb>`__)

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.6.1.svg
   :target: https://pypi.org/project/bezier/0.6.1/
   :alt: PyPI link to release 0.6.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.6.1
   :target: https://bezier.readthedocs.io/en/0.6.1/
   :alt: Documentation for release 0.6.1
