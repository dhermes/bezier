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
   CPython versions: 2.7, 3.6 and 3.7 on Linux, Mac OS X and Windows.
   In addition, the feature should have 100% line coverage.

#. **Documentation**: The feature must (should) be documented with
   helpful `doctest`_ examples wherever relevant.

.. _File an issue: https://github.com/dhermes/bezier/issues/new
.. _doctest: http://www.sphinx-doc.org/en/stable/ext/doctest.html

****************
Binary Extension
****************

Many low-level computations have alternate implementations in Fortran.
See the `Python Binary Extension`_ page in the documentation for a more
detailed description.

.. _Python Binary Extension: https://bezier.readthedocs.io/en/latest/python/binary-extension.html

Building / Compiling
====================

To compile the binary extension (with a Fortran compiler installed) run:

.. code-block:: console

   $ python setup.py build_ext --inplace
   $ # OR
   $ python setup.py build_ext --inplace --fcompiler=${FC}

By default the Fortran code will be compiled in "optimized" mode, which
may make debugging more difficult. To compile with debug symbols, without
optimizations that move code around, etc. use the ``DEBUG`` environment
variable:

.. code-block:: console

   $ DEBUG=True python setup.py build_ext --inplace

Using ``distutils`` and ``numpy.distutils`` to compile Fortran is not
"fully-supported" (i.e. the tooling is ad-hoc). As a result, there is a
decent amount of code in ``setup.py``, ``setup_helpers.py``,
``setup_helpers_osx.py`` and ``setup_helpers_windows.py`` to specify the build
process. To make sure these are working as expected, it's possible to
track **how** extensions are being installed. To actually make sure the
correct compiler commands are invoked, provide a filename as the
``BEZIER_JOURNAL`` environment variable and then the commands invoked will
be written there:

.. code-block:: console

   $ BEZIER_JOURNAL=$(pwd)/journal.txt python setup.py build_ext --inplace

The ``nox`` session ``check_journal`` uses this journaling option to verify
the commands used to compile the extensions in Linux on `CircleCI`_, in
Mac OS X on `Travis CI`_ and in Windows on `AppVeyor`_.

As the build complexity grows, it may make more sense to transition the steps
out of Python and into `CMake`_, `SCons`_ or another build tool.

.. _CMake: https://cmake.org
.. _SCons: http://scons.org

To explicitly disable the building of extensions, the ``BEZIER_NO_EXTENSIONS``
environment variable can be used:

.. code-block:: console

   $ BEZIER_NO_EXTENSIONS=True .../bin/python -m pip install .

This environment variable is actually used for the ``nox -s docs`` session
to emulate the `RTD`_ build environment (where no Fortran compiler is
present).

Dependencies
============

Currently, the ``src/bezier/quadpack`` `directory`_ has a subset of Fortran 77
subroutines from `QUADPACK`_. These are Public Domain, so they do not
conflict with the Apache 2.0 license (as far as we know). In addition it
contains another popular subroutine from NETLIB: ``d1mach`` (which the
QUADPACK subroutines depend on).

QUADPACK is used to perform numerical quadrature to compute the length
of a curve segment.

.. _directory: https://github.com/dhermes/bezier/tree/master/src/bezier/quadpack
.. _QUADPACK: https://en.wikipedia.org/wiki/QUADPACK

******************
Running Unit Tests
******************

We recommend using `Nox`_ to run unit tests:

.. code-block:: console

   $ nox -s "unit-2.7"
   $ nox -s "unit-3.6"
   $ nox -s "unit-3.7"
   $ nox -s "unit-pypy"
   $ nox -s  unit  # Run all versions

However, `pytest`_ can be used directly (though it won't
manage dependencies or build extensions):

.. code-block:: console

   $ PYTHONPATH=src/ python2.7 -m pytest tests/unit/
   $ PYTHONPATH=src/ python3.6 -m pytest tests/unit/
   $ PYTHONPATH=src/ python3.7 -m pytest tests/unit/
   $ PYTHONPATH=src/ pypy      -m pytest tests/unit/

.. _Nox: https://nox.readthedocs.io
.. _pytest: https://docs.pytest.org

Testing the Binary Extension
============================

When using ``nox``, the ``bezier`` package will automatically be installed
into a virtual environment and the binary extension will be built during
install.

However, if the tests are run directly from the source tree via

.. code-block:: console

   $ PYTHONPATH=src/ python -m pytest tests/unit/

some unit tests may be skipped. The unit tests that explicitly exercise the
binary extension will skip (rather than fail) if the extension isn't
compiled (with ``build_ext --inplace``) and present in the source tree.

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
   $ PYTHONPATH=src/ python -m pytest \
   >   --cov=bezier \
   >   --cov=tests.unit \
   >   tests/unit/ \
   >   tests/functional/test_segment_box.py

Slow Tests
==========

To run unit tests without tests that have been (explicitly)
marked slow, use the ``--ignore-slow`` flag:

.. code-block:: console

   $ nox -s "unit-2.7" -- --ignore-slow
   $ nox -s "unit-3.6" -- --ignore-slow
   $ nox -s "unit-3.7" -- --ignore-slow
   $ nox -s  unit      -- --ignore-slow

These slow tests have been identified via:

.. code-block:: console

   $ ...
   $ nox -s "unit-3.7" -- --durations=10

and then marked with ``pytest.mark.skipif``.

Slow Install
============

Installing NumPy with `PyPy`_ can take upwards of two minutes and
installing SciPy can take as much as seven minutes. This makes it
prohibitive to create a new environment for testing.

.. _PyPy: https://pypy.org/

In order to avoid this penalty, the ``WHEELHOUSE`` environment
variable can be used to instruct ``nox`` to install NumPy and SciPy
from locally built wheels when installing the ``pypy`` sessions.

To pre-build NumPy and SciPy wheels:

.. code-block:: console

   $ pypy -m virtualenv pypy-venv
   $ pypy-venv/bin/python -m pip wheel --wheel-dir=${WHEELHOUSE} numpy
   $ pypy-venv/bin/python -m pip install ${WHEELHOUSE}/numpy*.whl
   $ pypy-venv/bin/python -m pip wheel --wheel-dir=${WHEELHOUSE} scipy
   $ rm -fr pypy-venv/

Alternatively, wheels can be downloaded from `pypy-wheels`_, however
the SciPy wheel will still require ``libatlas-dev``, ``libblas-dev`` and
``liblapack-dev``.

The `Docker`_ image for the CircleCI test environment has already
pre-built these wheels and stored them in the ``/wheelhouse`` directory.
So, in the `CircleCI environment`_, the ``WHEELHOUSE`` environment
variable is set to ``/wheelhouse``.

.. _Docker: https://www.docker.com/
.. _pypy-wheels: https://antocuni.github.io/pypy-wheels/

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

   $ nox -s "functional-2.7"
   $ nox -s "functional-3.6"
   $ nox -s "functional-3.7"
   $ nox -s "functional-pypy"
   $ nox -s  functional  # Run all versions
   $ # OR
   $ PYTHONPATH=src/ python2.7 -m pytest tests/functional/
   $ PYTHONPATH=src/ python3.6 -m pytest tests/functional/
   $ PYTHONPATH=src/ python3.7 -m pytest tests/functional/
   $ PYTHONPATH=src/ pypy      -m pytest tests/functional/

.. _functional tests: https://github.com/dhermes/bezier/tree/master/tests/functional

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
``tests/functional/schema`` directory.

.. _curves.json: https://github.com/dhermes/bezier/blob/master/tests/functional/curves.json
.. _curve_intersections.json: https://github.com/dhermes/bezier/blob/master/tests/functional/curve_intersections.json
.. _surfaces.json: https://github.com/dhermes/bezier/blob/master/tests/functional/surfaces.json
.. _surface_intersections.json: https://github.com/dhermes/bezier/blob/master/tests/functional/surface_intersections.json
.. _JSON schema: http://json-schema.org/

************
Coding Style
************

Code is `PEP8`_ compliant and this is enforced with `flake8`_
and `Pylint`_.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _flake8: http://flake8.pycqa.org
.. _Pylint: https://www.pylint.org

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
.. _Napoleon: https://sphinxcontrib-napoleon.readthedocs.io
.. _sphinx-docstring-typing: https://pypi.org/project/sphinx-docstring-typing/
.. _type annotation: https://docs.python.org/3/library/typing.html

*************
Documentation
*************

The documentation is built with `Sphinx`_ and automatically
updated on `RTD`_ every time a commit is pushed to ``master``.

.. _Sphinx: http://www.sphinx-doc.org
.. _RTD: https://readthedocs.org/

To build the documentation locally:

.. code-block:: console

   $ nox -s docs
   $ # OR (from a Python 3.6 or later environment)
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
   $ # OR (from a Python 3.6 or later environment)
   $ PYTHONPATH=src/ sphinx-build -W \
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
   $ # OR (from a Python 3.6 or later environment)
   $ export MATPLOTLIBRC=docs/ GENERATE_IMAGES=True PYTHONPATH=src/
   $ sphinx-build -W \
   >   -b doctest \
   >   -d docs/build/doctrees \
   >   docs \
   >   docs/build/doctest
   $ python tests/functional/test_segment_box.py --save-plot
   $ python tests/functional/test_surface_locate.py --save-plot
   $ python tests/functional/make_curve_curve_images.py
   $ python tests/functional/make_surface_surface_images.py
   $ unset MATPLOTLIBRC GENERATE_IMAGES PYTHONPATH

**********************
Continuous Integration
**********************

Tests are run on `CircleCI`_ (Linux), `Travis CI`_ (Mac OS X) and
`AppVeyor`_ (Windows) after every commit. To see which tests are run, see
the `CircleCI config`_, the `Travis config`_ and the `AppVeyor config`_.

On CircleCI, a `Docker`_ image is used to provide fine-grained control over
the environment. There is a base `python-multi Dockerfile`_ that just has the
Python versions we test in. The image used in our CircleCI builds (from
`bezier Dockerfile`_) installs dependencies needed for testing (such as
``nox`` and NumPy).

On Travis CI, Matthew Brett's `multibuild`_ is used to install "official"
python.org CPython binaries for Mac OS X. Then tests are run in both 64-bit
mode (NumPy has `discontinued`_ 32-bit support).

On AppVeyor, all extensions are built and tested with both 32-bit and 64-bit
Python binaries.

.. _CircleCI: https://circleci.com/gh/dhermes/bezier
.. _Travis CI: https://travis-ci.org/dhermes/bezier
.. _AppVeyor: https://ci.appveyor.com/project/dhermes/bezier
.. _CircleCI config: https://github.com/dhermes/bezier/blob/master/.circleci/config.yml
.. _Travis config: https://github.com/dhermes/bezier/blob/master/.travis.yml
.. _AppVeyor config: https://github.com/dhermes/bezier/blob/master/.appveyor.yml
.. _python-multi Dockerfile: https://github.com/dhermes/python-multi/blob/master/src/Dockerfile
.. _bezier Dockerfile: https://github.com/dhermes/bezier/blob/master/scripts/docker/bezier.Dockerfile
.. _multibuild: https://github.com/matthew-brett/multibuild
.. _discontinued: https://github.com/numpy/numpy/issues/11625

****************************************
Release Process / Deploying New Versions
****************************************

New versions are pushed to `PyPI`_ manually after a ``git`` tag is
created. The process is manual (rather than automated) for several
reasons:

* The documentation and README (which acts as the landing page text on
  PyPI) will be updated with links scoped to the versioned tag (rather
  than ``master``). This update occurs via the ``doc_template_release.py``
  script.
* Several badges on the documentation landing page (``index.rst``) are
  irrelevant to a fixed version (such as the "latest" version of the
  package).
* The build badges in the README and the documentation will be
  changed to point to a fixed (and passing) build that has already
  completed (will be the build that occurred when the tag was pushed). If
  the builds pushed to PyPI automatically, a build would need to
  link to itself **while** being run.
* Wheels need be built for Linux, Mac OS X and Windows. This process
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
-  `Python 3.6`_
-  `Python 3.7`_

.. _Python 2.7: https://docs.python.org/2.7/
.. _Python 3.6: https://docs.python.org/3.6/
.. _Python 3.7: https://docs.python.org/3.7/

Supported versions can be found in the ``noxfile.py`` `config`_.

.. _config: https://github.com/dhermes/bezier/blob/master/noxfile.py

Versioning
==========

``bezier`` follows `semantic versioning`_.

.. _semantic versioning: http://semver.org/

It is currently in major version zero (``0.y.z``), which means that
anything may change at any time and the public API should not be
considered stable.

*********************
Environment Variables
*********************

This project uses environment variables for building the
``bezier._speedup`` binary extension:

- ``BEZIER_JOURNAL``: If set to a path on the filesystem, all compiler
  commands executed while building the binary extension will be logged to
  the journal file
- ``BEZIER_NO_EXTENSIONS``: If set, this will indicate that only the pure
  Python package should be built and installed (i.e. without the binary
  extension).
- ``BEZIER_WHEEL``: Indicates that the source is being built into a wheel.
  When this is true, some compiler flags (e.g. ``-march=native``) will be
  removed since those flags can produce machine instructions that are too
  specific to the host platform / architecture.
- ``DEBUG``: Indicates the binary extension should be built in debug mode.

for interacting with the system at import time:

- ``PATH``: On Windows, we add the ``bezier/extra-dll`` package directory to
  the path so that the ``libbezier.dll`` shared libary can be loaded at
  import time

and for running tests and interacting with Continuous Integration
services:

- ``WHEELHOUSE``: If set, this gives a path to prebuilt NumPy and SciPy wheels
  for PyPy.
- ``GENERATE_IMAGES``: Indicates to ``nox -s doctest`` that images should
  be generated during cleanup of each test case.
- ``APPVEYOR``: Indicates currently running on AppVeyor.
- ``CIRCLECI``: Indicates currently running on CircleCI.
- ``READTHEDOCS``: Indicates currently running on Read The Docs (RTD). This is
  used to tell Sphinx to use the RTD theme when **not** running on RTD.
- ``TRAVIS``: Indicates currently running on Travis.
- ``TRAVIS_BUILD_DIR``: Gives path to the Travis build directory. This is used
  to modify the command journal to make it deterministic (i.e. independent
  of the build directory).
- ``TRAVIS_OS_NAME``: Gives the current operating system on Travis. We check
  that it is ``osx``.
