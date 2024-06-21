Latest Release (``2024.6.20``)
==============================

|pypi| |docs|

``libbezzier`` Changes
----------------------

Packaging
~~~~~~~~~

-  Stop building with ``-static`` on Windows
   (`#311 <https://github.com/dhermes/bezier/pull/311>`__).
   This is entirely enabled by the amazing
   `delvewheel <https://github.com/adang1345/delvewheel>`__ project (non-static
   dependencies can now easily be packaged into a built wheel for Python).

Python Changes
--------------

Packaging
~~~~~~~~~

-  Added support for Python 3.12
   (`#315 <https://github.com/dhermes/bezier/pull/315>`__).
-  Dropped support for Python 3.8
   (`#310 <https://github.com/dhermes/bezier/pull/310>`__)
   and Python 3.9
   (`#315 <https://github.com/dhermes/bezier/pull/315>`__).

Documentation
-------------

-  Add ``DEVELOPMENT.rst`` section about configuring a shell and compilers for
   Windows development
   (`#311 <https://github.com/dhermes/bezier/pull/311>`__). This had been a
   hodge podge of local and remote (CI) development for the last 5+ years, so
   this a big milestone!

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/2024.6.20.svg
   :target: https://pypi.org/project/bezier/2024.6.20/
   :alt: PyPI link to release 2024.6.20
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=2024.6.20
   :target: https://bezier.readthedocs.io/en/2024.6.20/
   :alt: Documentation for release 2024.6.20
