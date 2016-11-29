###########
Development
###########


.. contents:: Here are some pointers for hacking on ``bezier``.

***************
Adding Features
***************

In order to add a feature to ``bezier``:

#. **Discuss**: `File an issue`_ to notify maintainers of the
   proposed changes (i.e. just sending a large PR with a finished
   feature may catch maintainers off guard)

#. **Add tests**: The feature must work fully on the following
   CPython versions: 2.7, 3.4 and 3.5 on both UNIX and Windows.
   In addition, the feature should have 100% line coverage.

#. **Documentation**: The feature must (should) be documented with
   helpful `doctest`_ examples wherever relevant.

.. _File an issue: https://github.com/dhermes/bezier/issues/new
.. _doctest: http://www.sphinx-doc.org/en/stable/ext/doctest.html

******************
Running Unit Tests
******************

We recommend using ``tox`` to run unit tests:

.. code-block:: console

   $ tox -e py27
   $ tox -e py34
   $ tox -e py35

However, `pytest`_ can be used directly (though it won't
manage dependencies):

.. _pytest: http://docs.pytest.org/en/stable/

.. code-block:: console

   $ python2.7 -m py.test tests/
   $ python3.4 -m py.test tests/
   $ python3.5 -m py.test tests/

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

   $ tox -e cover
   $ # OR
   $ python -m py.test \
   >  --cov=bezier \
   >  --cov=tests \
   >  tests/ \
   >  functional_tests/

Slow Tests
==========

To run unit tests without tests that have been (explicitly)
marked slow, use the ``--ignore-slow`` flag:

.. code-block:: console

   $ python -m py.test tests/ --ignore-slow

These slow tests have been identified via:

.. code-block:: console

   $ python -m py.test tests/ --durations=10

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

   $ tox -e functional
   $ # OR
   $ python -m py.test functional_tests/

.. _functional tests: https://github.com/dhermes/bezier/tree/master/functional_tests

For example, the following curve-curve intersection is a
functional test case:

.. image:: https://github.com/dhermes/bezier/blob/master/docs/images/test_curves11_and_26.png

and there is a `curve-curve doc`_ which captures most of the cases in the
functional tests.

.. _curve-curve doc: http://bezier.readthedocs.io/en/latest/curve-curve-intersection.html

A surface-surface intersection functional test case:

.. image:: https://github.com/dhermes/bezier/blob/master/docs/images/test_surfaces1Q_and_2Q.png

and a segment-box functional test case:

.. image:: https://github.com/dhermes/bezier/blob/master/docs/images/test_goes_through_box.png

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

   $ tox -e lint

A few extensions and overrides have been specified in the `pylintrc`_
configuration for ``bezier``.

.. _pylintrc: https://github.com/dhermes/bezier/blob/master/pylintrc

Docstring Style
===============

We require docstrings on all public objects and enforce this with
our ``lint`` checks. The docstrings mostly follow `PEP257`_
and are written in the `Google style`_, e.g.

.. code-block::

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
.. _Google style: http://google.github.io/styleguide/pyguide.html#Comments__body
.. _Napoleon: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/
.. _sphinx-docstring-typing: https://pypi.python.org/pypi/sphinx-docstring-typing
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

   $ tox -e docs
   $ # OR
   $ ./scripts/build_docs.sh

Documentation Snippets
======================

A large effort is made to provide useful snippets in documentation.
To make sure these snippets are valid (and remain valid over
time), `doctest`_ is used to check that the interpreter output
in the snippets are valid.

To run the documentation tests:

.. code-block:: console

   $ tox -e doctest
   $ # OR
   $ sphinx-build -W \
   >   -b doctest \
   >   -d docs/build/doctrees \
   >   docs \
   >   build/doctest

**********************
Continuous Integration
**********************

And

**********************
Deploying New Versions
**********************

New versions are deployed to `PyPI`_ automatically every time
a new tag is pushed. To allow `twine`_ to authenticate (which
is needed to upload) the ``TWINE_USERNAME`` and ``TWINE_PASSWORD``
environment variables are set in the `CircleCI environment`_.

.. _PyPI: https://pypi.python.org/pypi/bezier
.. _twine: https://packaging.python.org/distributing/

Supported Python Versions
=========================

``bezier`` explicitly supports:

-  `Python 2.7`_
-  `Python 3.4`_
-  `Python 3.5`_

.. _Python 2.7: https://docs.python.org/2.7/
.. _Python 3.4: https://docs.python.org/3.4/
.. _Python 3.5: https://docs.python.org/3.5/

Supported versions can be found in the ``tox.ini`` `config`_.

.. _config: https://github.com/dhermes/bezier/blob/master/tox.ini

Versioning
==========

``bezier`` follows `semantic versioning`_.

.. _semantic versioning: http://semver.org/

It is currently in major version zero (``0.y.z``), which means that
anything may change at any time and the public API should not be
considered stable.
