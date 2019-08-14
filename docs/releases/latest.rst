Latest Release (``0.10.0``)
===========================

|pypi| |docs|

Python Changes
--------------

Breaking Changes
~~~~~~~~~~~~~~~~

-  Support for Python 2.7 has been dropped. With the impending `EOL`_ of Python
   2.7 on January 1, 2020 many of the ``bezier`` dependencies such as
   ``numpy``, ``scipy`` and ``pytest`` have dropped support for Python 2.7 in
   their latest releases. Some changes related to this include:

   -  Removing support for Python 2.7
      (`3eaa5aa <https://github.com/dhermes/bezier/commit/3eaa5aaa670d167b2c1340d3d531d5438eaf62cd>`__).
   -  Updating all PyPy usage to be Python 3 only
      (`1e3037f <https://github.com/dhermes/bezier/commit/1e3037fce5acdcfa194cac481ee06ef6bcc329e5>`__).
   -  Removing ``_bool_patch.h``
      (`4ccc559 <https://github.com/dhermes/bezier/commit/4ccc559e6928f78556c1201f45a2ad7b3b40d7a5>`__).

.. _EOL: https://pythonclock.org/

Build
~~~~~

-  Integrating ``black`` `code formatter`_ into the ``nox -s lint`` session
   to ensure consistent formatting
   (`e659532 <https://github.com/dhermes/bezier/commit/e659532747d0433bf3a91198a7baf172ed36f069>`__).
-  Fully automating the building of wheels in the ``bezier-wheels`` `project`_
   (`recent commits`_). Built wheels are uploaded to a Google Cloud Storage
   bucket.
-  Using the same set of optimal flags for Fortran 77 code (i.e. ``.f`` files)
   that are used for Fortran 90 code (i.e. ``.f90`` files)
   (`e7eb56e <https://github.com/dhermes/bezier/commit/e7eb56e723f13d43f6eae855e6556b4ccbc1edd9>`__).
-  Unify ``requirements.txt`` files and add notes about why each dependency is
   required
   (`230814d <https://github.com/dhermes/bezier/commit/230814d67e24f42f967a652ff7e8d81ee2176954>`__,
   `1ae147f <https://github.com/dhermes/bezier/commit/1ae147f81e7a01ba672806a8fd56de25ba2bdcdb>`__,
   `e710ee6 <https://github.com/dhermes/bezier/commit/e710ee6968438cb2462ec8bea8af407159a63925>`__).
-  Changing ``imp`` usage to ``importlib`` due to deprecation of the former
   (`9231d92 <https://github.com/dhermes/bezier/commit/9231d92b420df1ed97ae2b159bd0aedf0c1ff888>`__).
   Fixed `#137 <https://github.com/dhermes/bezier/issues/137>`__.
-  Ditching the ``--check-archs`` flag in the macOS script for building wheels
   since we can no longer support 32-bit on macOS due to NumPy
   (`37be384 <https://github.com/dhermes/bezier/commit/37be3845750ff0fe9f200f87a8427b05639c3a61>`__).
-  Improved dev experience with Docker image used on CircleCI by adding a
   ``.dockerignore`` file for faster builds, suggesting ``--workdir`` and
   flag during local dev, setting ``WHEELHOUSE`` environment variable directly
   in the container (rather than in the CircleCI settings) and allowing
   "default" locations for pre-built wheels at ``/wheelhouse`` and
   ``${HOME}/wheelhouse``
   (`08be336 <https://github.com/dhermes/bezier/commit/08be336efac467beeb7055cfc80996b97482456a>`__,
   `26acc38 <https://github.com/dhermes/bezier/commit/26acc384d857cf9f5ddd8260ef50b7bcffeeb133>`__,
   `7634779 <https://github.com/dhermes/bezier/commit/763477958c73a4eb6ce0f89b6b37887c66c10706>`__,
   `f9a8fcf <https://github.com/dhermes/bezier/commit/f9a8fcf275b244d962fae1e93b223af0c78285cc>`__).

.. _recent commits: https://github.com/dhermes/bezier-wheels/compare/ee008511d5ff2736dfb44f770552e7553b00e8f0...424453f50fbb8f240ca60280b637a278f6e9ad4a
.. _code formatter: https://black.readthedocs.io

Miscellany
~~~~~~~~~~

-  Make some functional test cases more lenient so that they pass on 32-bit
   CentOS 5, which is used for ``manylinux``
   (`e7eb56e <https://github.com/dhermes/bezier/commit/e7eb56e723f13d43f6eae855e6556b4ccbc1edd9>`__).
   This was part of a large effort to fully automate the building of wheels in
   the ``bezier-wheels`` `project`_.
-  Replacing ``pypy`` with ``pypy3`` in testing as the only non-CPython
   "unofficially supported" runtime. (This is part of the drop in support for
   Python 2.7.) Unfortunately the currently (as of August 2019) released
   versions of ``pypy3`` are not currently working with ``numpy >= 1.16``
   (see `numpy/numpy#12740 <https://github.com/numpy/numpy/issues/12740>`__)
   so the ``numpy == 1.15.4`` version is a pinned dependency.

   -  Specifying the NumPy version in ``setup.py`` based on
      ``implementation_name``
      (`7e9046d <https://github.com/dhermes/bezier/commit/7e9046dc9dbe6f448238141221c5a7dff497d8d4>`__).
   -  Add ``_pypy_speedup.c`` built with Cython 0.29.11 because the latest
      Cython (0.29.13 as of August 2019) corresponds to the versions of NumPy
      that are incompatible with PyPy
      (`7813e41 <https://github.com/dhermes/bezier/commit/7813e41f7666fa36fbb4a7daf0aa45c2d2bee87f>`__).
   -  Pinning to ``numpy==1.15.4`` and ``scipy==1.2.0`` in wheelhouse for
      pre-built Docker container
      (`7634779 <https://github.com/dhermes/bezier/commit/763477958c73a4eb6ce0f89b6b37887c66c10706>`__).

-  Added ``nox -s validate_functional_test_cases`` session to ensure that
   functional test cases always adhere to the JSON schema.

   -  Added ``nox`` session and fixed some schema file bugs
      (`618653a <https://github.com/dhermes/bezier/commit/618653a0888cc5e91a5fb1959cf5e04f61e5c1cf>`__).
   -  Fixed curve intersections that did not adhere to schema
      (`35b158a <https://github.com/dhermes/bezier/commit/35b158a9ad4f8c0ed1d4a3cd07a8c157f33b0639>`__).
   -  Transposed ``nodes`` schema for curved polygon
      (`8c3ca89 <https://github.com/dhermes/bezier/commit/8c3ca895512a60c2fe82d8a24ab328244e3abb3f>`__).
   -  Enable the ``nox`` session to run in CircleCI
      (`5a0a343 <https://github.com/dhermes/bezier/commit/5a0a343728ac52933b1aadd3c483fb439f2e043a>`__).

-  Updated ``@slow`` marker for ``pytest`` because it used a deprecated API
   (`46f8b57 <https://github.com/dhermes/bezier/commit/46f8b57c8b34484236ce1bc9aa9f5ea5fc77c5df>`__).

.. _project: https://github.com/dhermes/bezier-wheels

Documentation
-------------

- Changing all references to Mac OS X to macOS
  (`b10c2fc <https://github.com/dhermes/bezier/commit/b10c2fc1af424e862143ac40d01f7baa65fc8af0>`__,
  `c1c2c6b <https://github.com/dhermes/bezier/commit/c1c2c6b767c40c2eb070ae599a110ecc9fb3e793>`__,
  `131d17b <https://github.com/dhermes/bezier/commit/131d17be3db5546deebff953378252b12b426534>`__).
  As of 10.12, the operating system has changed its name.
- Splitting up ``algorithms/helpers``. The pre-amble has been moved into the
  ``algorithms`` landing page and the geometric and algebraic helpers have been
  moved into separate docs.
  (`889c913 <https://github.com/dhermes/bezier/commit/889c913436b6d01533d8eb1830717620cea725ef>`__).

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.10.0.svg
   :target: https://pypi.org/project/bezier/0.10.0/
   :alt: PyPI link to release 0.10.0
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.10.0
   :target: https://bezier.readthedocs.io/en/0.10.0/
   :alt: Documentation for release 0.10.0
