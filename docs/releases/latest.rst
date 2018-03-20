Latest Release (``0.8.1.dev1``)
===============================

|pypi| |docs|

ABI Changes
-----------

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removing getters and setters for parameters used during curve-curve
   intersection
   (`2fda3ae <https://github.com/dhermes/bezier/commit/2fda3aed2818849363c425e3fce70b4bafe7e9ef>`__):

   -  ``curve_intersection.h::set_max_candidates``
   -  ``curve_intersection.h::get_max_candidates``

Python Changes
--------------

Non-Public API
~~~~~~~~~~~~~~

-  Removing getters and setters for parameters used during curve-curve
   intersection
   (`2fda3ae <https://github.com/dhermes/bezier/commit/2fda3aed2818849363c425e3fce70b4bafe7e9ef>`__):

   -  ``bezier._geometric_intersection.set_max_candidates()``
   -  ``bezier._geometric_intersection.get_max_candidates()``

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.8.1.svg
   :target: https://pypi.org/project/bezier/0.8.1/
   :alt: PyPI link to release 0.8.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.8.1
   :target: https://bezier.readthedocs.io/en/0.8.1/
   :alt: Documentation for release 0.8.1
