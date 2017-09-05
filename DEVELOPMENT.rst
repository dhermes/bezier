###########
Development
###########


.. contents:: Here are some pointers for hacking on ``bezier``.

***************
Adding Features
***************

In order to add a feature to ``bezier``:

#. **Discuss**: `File an issue`_ to notify maintainer(s) of the
   proposed changes (i.e. just sending a large PR with a finished
   feature may catch maintainer(s) off guard).

#. **Add tests**: The feature must work fully on the following
   CPython versions: 2.7, 3.5 and 3.6 on both UNIX and Windows.
   In addition, the feature should have 100% line coverage.

#. **Documentation**: The feature must (should) be documented with
   helpful `doctest`_ examples wherever relevant.

.. _File an issue: https://github.com/dhermes/bezier/issues/new
.. _doctest: http://www.sphinx-doc.org/en/stable/ext/doctest.html

******************
Running Unit Tests
******************

We recommend using ``nox`` (`nox-automation`_) to run unit tests:

.. code-block:: console

   $ nox -s "unit_tests(python_version='2.7')"
   $ nox -s "unit_tests(python_version='3.5')"
   $ nox -s "unit_tests(python_version='3.6')"
   $ nox -s "unit_tests(python_version='pypy')"
   $ nox -s unit_tests  # Run all versions

However, `pytest`_ can be used directly (though it won't
manage dependencies or build extensions):

.. code-block:: console

   $ PYTHONPATH=src/ python2.7 -m pytest tests/
   $ PYTHONPATH=src/ python3.5 -m pytest tests/
   $ PYTHONPATH=src/ python3.6 -m pytest tests/
   $ PYTHONPATH=src/ MATPLOTLIBRC=test/ pypy -m pytest tests/

.. _nox-automation: https://nox.readthedocs.io/en/stable/
.. _pytest: https://docs.pytest.org/en/stable/

Native Code Extensions
======================

Many low-level computations have alternate implementations in Fortran.
When using ``nox``, the ``bezier`` package will automatically be installed
into a virtual environment and the native extensions will be built during
install.

However, if you run the tests directly from the source tree via

.. code-block:: console

   $ PYTHONPATH=src/ python -m pytest tests/

some unit tests may be skipped. The unit tests for the Fortran
implementations will skip (rather than fail) if the extensions aren't
compiled and present in the source tree. To compile the native extensions,
make sure you have a valid Fortran compiler and run

.. code-block:: console

   $ python setup.py build_ext --inplace
   $ # OR
   $ python setup.py build_ext --inplace --fcompiler=${FC}

Test Coverage
=============

``bezier`` has 100% `line coverage`_. The coverage is checked
on every build and uploaded to `coveralls.io`_ via the
``COVERALLS_REPO_TOKEN`` environment variable set in
the `CircleCI environment`_.

.. _line coverage: https://coveralls.io/github/dhermes/bezier
.. _coveralls.io: https://coveralls.io/
.. _CircleCI environment: https://circleci.com/gh/dhermes/bezier/edit#env-vars

To run the coverage report locally:

.. code-block:: console

   $ nox -s cover
   $ # OR
   $ PYTHONPATH=src/:functional_tests/ python -m pytest \
   >   --cov=bezier \
   >   --cov=tests \
   >   tests/ \
   >   functional_tests/test_segment_box.py

Slow Tests
==========

To run unit tests without tests that have been (explicitly)
marked slow, use the ``--ignore-slow`` flag:

.. code-block:: console

   $ nox -s "unit_tests(python_version='2.7')" -- --ignore-slow
   $ nox -s "unit_tests(python_version='3.5')" -- --ignore-slow
   $ nox -s "unit_tests(python_version='3.6')" -- --ignore-slow
   $ nox -s  unit_tests                        -- --ignore-slow

These slow tests have been identified via:

.. code-block:: console

   $ ...
   $ nox -s "unit_tests(python_version='3.6')" -- --durations=10

and then marked with ``pytest.mark.skipif``.

****************
Functional Tests
****************

Line coverage and unit tests are not entirely sufficient to
test **numerical software**. As a result, there is a fairly
large collection of `functional tests`_ for ``bezier``.

These give a broad sampling of curve-curve intersection,
surface-surface intersection and segment-box intersection problems to
check both the accuracy (i.e. detecting all intersections) and the
precision of the detected intersections.

To run the functional tests:

.. code-block:: console

   $ nox -s "functional(python_version='2.7')"
   $ nox -s "functional(python_version='3.5')"
   $ nox -s "functional(python_version='3.6')"
   $ nox -s "functional(python_version='pypy')"
   $ nox -s  functional  # Run all versions
   $ # OR
   $ export PYTHONPATH=src/:functional_tests/
   $ python2.7 -m pytest functional_tests/
   $ python3.5 -m pytest functional_tests/
   $ python3.6 -m pytest functional_tests/
   $ MATPLOTLIBRC=test/ pypy -m pytest functional_tests/
   $ unset PYTHONPATH

.. _functional tests: https://github.com/dhermes/bezier/tree/master/functional_tests

For example, the following curve-curve intersection is a
functional test case:

.. image:: https://raw.githubusercontent.com/dhermes/bezier/master/docs/images/curves11_and_26.png
   :align: center

and there is a `Curve-Curve Intersection`_ document which captures many of
the cases in the functional tests.

.. _Curve-Curve Intersection: https://bezier.readthedocs.io/en/latest/curve-curve-intersection.html

A surface-surface intersection functional test case:

.. image:: https://raw.githubusercontent.com/dhermes/bezier/master/docs/images/surfaces1Q_and_2Q.png
   :align: center

a segment-box functional test case:

.. image:: https://raw.githubusercontent.com/dhermes/bezier/master/docs/images/test_goes_through_box08.png
   :align: center

and a "locate point on surface" functional test case:

.. image:: https://raw.githubusercontent.com/dhermes/bezier/master/docs/images/test_surface3_and_point1.png
   :align: center

Functional Test Data
====================

The curve-curve and surface-surface intersection test cases are stored in
JSON files:

* `curves.json`_
* `curve_intersections.json`_
* `surfaces.json`_
* `surface_intersections.json`_

This way, the test cases are programming language agnostic and can be
repurposed. The `JSON schema`_ for these files are stored in the
``functional_tests/schema`` directory.

.. _curves.json: https://github.com/dhermes/bezier/blob/master/functional_tests/curves.json
.. _curve_intersections.json: https://github.com/dhermes/bezier/blob/master/functional_tests/curve_intersections.json
.. _surfaces.json: https://github.com/dhermes/bezier/blob/master/functional_tests/surfaces.json
.. _surface_intersections.json: https://github.com/dhermes/bezier/blob/master/functional_tests/surface_intersections.json
.. _JSON schema: http://json-schema.org/

************
Coding Style
************

Code is `PEP8`_ compliant and this is enforced with `flake8`_
and `pylint`_.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _flake8: http://flake8.pycqa.org/en/stable/
.. _pylint: https://www.pylint.org/

To check compliance:

.. code-block:: console

   $ nox -s lint

A few extensions and overrides have been specified in the `pylintrc`_
configuration for ``bezier``.

.. _pylintrc: https://github.com/dhermes/bezier/blob/master/pylintrc

Docstring Style
===============

We require docstrings on all public objects and enforce this with
our ``lint`` checks. The docstrings mostly follow `PEP257`_
and are written in the `Google style`_, e.g.

.. code-block:: rest

   Args:
       path (str): The path of the file to wrap
       field_storage (FileStorage): The :class:`FileStorage` instance to wrap
       temporary (bool): Whether or not to delete the file when the File
          instance is destructed

   Returns:
       BufferedFileStorage: A buffered writable file descriptor

In order to support these in Sphinx, we use the `Napoleon`_ extension.
In addition, the `sphinx-docstring-typing`_ Sphinx extension is used to
allow for `type annotation`_ for arguments and result (introduced in
Python 3.5).

.. _PEP257: https://www.python.org/dev/peps/pep-0257/
.. _Google style: https://google.github.io/styleguide/pyguide.html#Comments__body
.. _Napoleon: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/
.. _sphinx-docstring-typing: https://pypi.org/project/sphinx-docstring-typing/
.. _type annotation: https://docs.python.org/3/library/typing.html

*************
Documentation
*************

The documentation is built with `Sphinx`_ and automatically
updated on `RTD`_ every time a commit is pushed to ``master``.

.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _RTD: https://readthedocs.org/

To build the documentation locally:

.. code-block:: console

   $ nox -s docs
   $ # OR (from a Python 3.5 or later environment)
   $ PYTHONPATH=src/ ./scripts/build_docs.sh

Documentation Snippets
======================

A large effort is made to provide useful snippets in documentation.
To make sure these snippets are valid (and remain valid over
time), `doctest`_ is used to check that the interpreter output
in the snippets are valid.

To run the documentation tests:

.. code-block:: console

   $ nox -s doctest
   $ # OR (from a Python 3.5 or later environment)
   $ PYTHONPATH=src/ NO_IMAGES=True sphinx-build -W \
   >   -b doctest \
   >   -d docs/build/doctrees \
   >   docs \
   >   docs/build/doctest

Documentation Images
====================

Many images are included to illustrate the curves / surfaces / etc.
under consideration and to display the result of the operation
being described. To keep these images up-to-date with the doctest
snippets, the images are created as doctest cleanup.

In addition, the images in the `Curve-Curve Intersection`_ document and
this document are generated as part of the functional tests.

To regenerate all the images:

.. code-block:: console

   $ nox -s docs_images
   $ # OR (from a Python 3.5 or later environment)
   $ export MATPLOTLIBRC=docs/ PYTHONPATH=src/
   $ sphinx-build -W \
   >   -b doctest \
   >   -d docs/build/doctrees \
   >   docs \
   >   docs/build/doctest
   $ python functional_tests/test_segment_box.py --save-plot
   $ python functional_tests/test_surface_locate.py --save-plot
   $ python functional_tests/make_curve_curve_images.py
   $ python functional_tests/make_surface_surface_images.py
   $ unset MATPLOTLIBRC PYTHONPATH

**********************
Continuous Integration
**********************

Tests are run on `CircleCI`_ and `AppVeyor`_ after every commit. To see
which tests are run, see the `CircleCI config`_ and the
`AppVeyor config`_.

.. _CircleCI: https://circleci.com/gh/dhermes/bezier
.. _AppVeyor: https://ci.appveyor.com/project/dhermes/bezier
.. _CircleCI config: https://github.com/dhermes/bezier/blob/master/.circleci/config.yml
.. _AppVeyor config: https://github.com/dhermes/bezier/blob/master/.appveyor.yml

**********************
Deploying New Versions
**********************

New versions are pushed to `PyPI`_ manually after a ``git`` tag is
created. The process is manual (rather than automated) for several
reasons:

* The documentation and README (which acts as the landing page text on
  PyPI) will be updated with links scoped to the versioned tag (rather
  than ``master``).
* Several badges on the documentation landing page (``index.rst``) are
  irrelevant to a fixed version (such as the "latest" version of the
  package).
* The build badges in the README and the documentation will be
  changed to point to a fixed (and passing) build that has already
  completed (will be the build that occurred when the tag was pushed). If
  the builds pushed to PyPI automatically, a build would need to
  link to itself **while** being run.
* Wheels need be built for Linux, Windows and OS X. This process
  is **becoming** better, but is still scattered across many
  different build systems. Each wheel will be pushed directly to
  PyPI via `twine`_.
* The release will be manually pushed to `TestPyPI`_ so the landing
  page can be visually inspected and the package can be installed
  from TestPyPI rather than from a local file.

.. _PyPI: https://pypi.org/project/bezier/
.. _twine: https://packaging.python.org/distributing/
.. _TestPyPI: https://packaging.python.org/guides/using-testpypi/

Supported Python Versions
=========================

``bezier`` explicitly supports:

-  `Python 2.7`_
-  `Python 3.5`_
-  `Python 3.6`_

.. _Python 2.7: https://docs.python.org/2.7/
.. _Python 3.5: https://docs.python.org/3.5/
.. _Python 3.6: https://docs.python.org/3.6/

Supported versions can be found in the ``nox.py`` `config`_.

.. _config: https://github.com/dhermes/bezier/blob/master/nox.py

Versioning
==========

``bezier`` follows `semantic versioning`_.

.. _semantic versioning: http://semver.org/

It is currently in major version zero (``0.y.z``), which means that
anything may change at any time and the public API should not be
considered stable.
