Latest Release (``0.9.1.dev1``)
===============================

|pypi| |docs|

Documentation
-------------

- Changing all references to Mac OS X to macOS
  (`b10c2fc <https://github.com/dhermes/bezier/commit/b10c2fc1af424e862143ac40d01f7baa65fc8af0>`__,
  `c1c2c6b <https://github.com/dhermes/bezier/commit/c1c2c6b767c40c2eb070ae599a110ecc9fb3e793>`__,
  `131d17b <https://github.com/dhermes/bezier/commit/131d17be3db5546deebff953378252b12b426534>`__).
  As of 10.12, the operating system has changed its name.

Python Changes
--------------

Build
~~~~~

- Fully automating the building of wheels in the ``bezier-wheels`` `project`_
  (`recent commits`_). Built wheels are uploaded to a Google Cloud Storage
  bucket.
- Using the same set of optimal flags for Fortran 77 code (i.e. ``.f`` files)
  that are used for Fortran 90 code (i.e. ``.f90`` files)
  (`e7eb56e <https://github.com/dhermes/bezier/commit/e7eb56e723f13d43f6eae855e6556b4ccbc1edd9>`__).
- Unify ``requirements.txt`` files and add notes about why each dependency is
  required
  (`230814d <https://github.com/dhermes/bezier/commit/230814d67e24f42f967a652ff7e8d81ee2176954>`__,
  `1ae147f <https://github.com/dhermes/bezier/commit/1ae147f81e7a01ba672806a8fd56de25ba2bdcdb>`__,
  `e710ee6 <https://github.com/dhermes/bezier/commit/e710ee6968438cb2462ec8bea8af407159a63925>`__).
- Changing ``imp`` usage to ``importlib`` due to deprecation of the former
  (`9231d92 <https://github.com/dhermes/bezier/commit/9231d92b420df1ed97ae2b159bd0aedf0c1ff888>`__).
  Fixed `#137 <https://github.com/dhermes/bezier/issues/137>`__.
- Ditching the ``--check-archs`` flag in the macOS script for building wheels
  since we can no longer support 32-bit on macOS due to NumPy
  (`37be384 <https://github.com/dhermes/bezier/commit/37be3845750ff0fe9f200f87a8427b05639c3a61>`__).

.. _recent commits: https://github.com/dhermes/bezier-wheels/compare/ee008511d5ff2736dfb44f770552e7553b00e8f0...424453f50fbb8f240ca60280b637a278f6e9ad4a

Miscellany
~~~~~~~~~~

- Make some functional test cases more lenient so that they pass on 32-bit
  CentOS 5, which is used for ``manylinux``
  (`e7eb56e <https://github.com/dhermes/bezier/commit/e7eb56e723f13d43f6eae855e6556b4ccbc1edd9>`__).
  This was part of a large effort to fully automate the building of wheels in
  the ``bezier-wheels`` `project`_.

.. _project: https://github.com/dhermes/bezier-wheels

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.9.1.svg
   :target: https://pypi.org/project/bezier/0.9.1/
   :alt: PyPI link to release 0.9.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.9.1
   :target: https://bezier.readthedocs.io/en/0.9.1/
   :alt: Documentation for release 0.9.1
